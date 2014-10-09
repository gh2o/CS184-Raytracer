#include <getopt.h>
#include <png.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <locale>
#include <memory>
#include <map>
#include <vector>
#include <array>
#include <limits>
#include <Eigen/Core>
#include <Eigen/Geometry>

using namespace Eigen;

typedef Transform<double,3,Affine> Transform4d;
typedef Array<double,3,1> Color3d;

class ParseException : public std::runtime_error {
public:
	ParseException(std::string msg) :
		ParseException(msg, -1) {}
	ParseException(std::string msg, int lineno) :
		runtime_error(buildMessage(msg, lineno)) {}
	static void showWarning(std::string msg) {
		showWarning(msg, -1);
	}
	static void showWarning(std::string msg, int lineno) {
		std::cerr << "Warning: " << buildMessage(msg, lineno) << std::endl;
	}
private:
	static std::string buildMessage(std::string msg, int lineno) {
		if (lineno > 0) {
			std::ostringstream s;
			s << "line " << lineno << ": " << msg;
			return s.str();
		} else {
			return msg;
		}
	}
};

class MathException : public std::runtime_error {
public:
	using runtime_error::runtime_error;
};

class WriteException : public std::runtime_error {
public:
	using runtime_error::runtime_error;
};

class Util {
public:
	static Vector4d cross(const Vector4d& a, const Vector4d& b) {
		Vector4d c = Vector4d::Zero();
		c.head<3>() = a.head<3>().cross(b.head<3>());
		return c;
	}
	/* NOTE: vec?dFrom?d functions should only be used with
	 *       directional (non-homogeneous) vectors. */
	static Vector4d vec4dFrom3d(const Vector3d& f) {
		Vector4d t = Vector4d::Zero();
		t.head<3>() = f;
		return t;
	}
	static Vector3d vec3dFrom4d(const Vector4d& f) {
		return f.head<3>();
	}
};

class Ray {
public:
	Vector4d origin_;
	Vector4d direction_;
	Ray(Vector4d inputOrigin, Vector4d inputDirection) {
		if (inputDirection.isZero())
			throw MathException("ray has no direction");
		if (inputOrigin(3) == 0)
			throw MathException("ray origin is a direction vector");
		if (inputDirection(3) != 0)
			throw MathException("ray direction is a point vector");
		origin_ = inputOrigin;
		direction_ = inputDirection.normalized();
	}
};

class Material {
public:
	Color3d ambientColor_;
	Color3d diffuseColor_;
	Color3d specularColor_;
	Color3d reflectiveColor_;
	double specularCoefficient_;
};

class Transformable {
public:
	Transform4d transform_;
};

class Camera : public Transformable {
public:
	Vector4d eyePoint_;
	Vector4d lowerLeftPoint_;
	Vector4d lowerRightPoint_;
	Vector4d upperLeftPoint_;
	Vector4d upperRightPoint_;
};

class Light : public Transformable {
public:
	Color3d color_;
};

class PointLight : public Light {
public:
	enum Falloff {
		FALLOFF_NONE,
		FALLOFF_LINEAR,
		FALLOFF_QUADRATIC,
	};
	Vector4d point_;
	Falloff falloff_;
};

class DirectionalLight : public Light {
public:
	Vector4d direction_;
};

class AmbientLight : public Light {
};

class Geometry : public Transformable {
public:
	Material material_;

	bool calculateIntersectionNormal(Ray inputRay, Vector4d& intersectionPt, Vector4d& normalDirection) {
		Matrix4d fwdTransform = transform_.matrix();
		Matrix4d invTransform = transform_.inverse().matrix();
		Ray transformedRay(invTransform * inputRay.origin_, invTransform * inputRay.direction_);
		Vector4d transformedIntersectionPt, transformedNormalDirection;
		bool hasIntersection = calculateIntNormInObjSpace(transformedRay, transformedIntersectionPt, 
			transformedNormalDirection);
		if (!hasIntersection) {
			return false;
		}
		intersectionPt = fwdTransform * transformedIntersectionPt;
		normalDirection = invTransform.transpose() * transformedNormalDirection;
		if (fwdTransform.determinant() < 0) {
			normalDirection = -normalDirection;
		}
		return true;
	}

	virtual bool calculateIntNormInObjSpace(Ray inputRay, Vector4d& intersectionPt, Vector4d& normalDirection) = 0;
};

class Sphere : public Geometry {
public:
	Vector4d center_;
	float radius_;

	bool calculateIntNormInObjSpace(Ray inputRay, Vector4d& intersectionPt, Vector4d& normalDirection) {
		Vector4d originCenterDiff = inputRay.origin_ - center_;
		double a = inputRay.direction_.dot(inputRay.direction_);
		double b = 2 * inputRay.direction_.dot(originCenterDiff);
		double c = originCenterDiff.dot(originCenterDiff) - radius_*radius_;
		double discriminant = b*b - 4*a*c;
		if (discriminant < 0) {
			return false;
		}
		double leftSide = (-b) / (2*a);
		double rightSide = sqrt(discriminant) / (2*a);
		double resultLower = leftSide - rightSide;
		double resultUpper = leftSide + rightSide;
		double result;
		if (resultUpper < 0) {
			return false;
		}
		result = resultLower > 0 ? resultLower : resultUpper;
		intersectionPt = inputRay.origin_ + result * inputRay.direction_;
		normalDirection = intersectionPt - center_;
		return true;
	}
};


class Mesh : public Geometry {
public:
	struct Face { std::array<Vector4d,3> points_, normals_; };
	std::vector<Face> faces_;
public:
	void addTriangle(const std::array<Vector4d,3>& points) {
		Vector4d a = points[1] - points[0];
		Vector4d b = points[2] - points[0];
		Vector4d n = Util::cross(a, b).normalized();
		Face face;
		face.points_ = points;
		face.normals_ = {{ n, n, n }};
		faces_.push_back(face);
	}
	bool calculateIntNormInObjSpace(Ray inputRay, Vector4d& intersectionPt, Vector4d& normalDirection) {
		bool closestExists = false;
		double closestDist = std::numeric_limits<double>::infinity();
		for (Face& face : faces_) {
			// base vectors for matrix
			Vector4d va = face.points_[1] - face.points_[0];
			Vector4d vb = face.points_[2] - face.points_[0];
			Vector4d dir = inputRay.direction_;
			Matrix4d m;
			m << va, vb, -dir, Vector4d(0,0,0,1);
			if (m.determinant() == 0)
				continue;
			// solve matrix for va and vb coeffs, direction scale t
			Vector4d sol = m.householderQr().solve(inputRay.origin_ - face.points_[0]);
			double a = sol(0);
			double b = sol(1);
			double t = sol(2);
			// validate constraints
			if (a < 0 || b < 0 || a + b > 1)
				continue;
			if (t < 0)
				continue;
			// check if closest
			double dist = t * dir.norm();
			if (dist >= closestDist)
				continue;
			// update closest
			closestExists = true;
			closestDist = dist;
			intersectionPt = face.points_[0] + a * va + b * vb;
			normalDirection =
				(1.0 - a - b) * face.normals_[0] +
				a * face.normals_[1] +
				b * face.normals_[2];
		}
		return closestExists;
	}
};

class GlobalScene {
public:
	typedef Array<Vector3d, Dynamic, Dynamic, RowMajor> RasterImage;
public:
	GlobalScene() :
		hasCamera_(false) {}
	void renderScene(RasterImage& output) {
		for (int r = 0; r < output.rows(); r++) {
			for (int c = 0; c < output.cols(); c++) {
				double row = r / output.rows();
				double col = c / output.cols();
				Vector4d pointOnImagePlane = col * (row * camera_.lowerLeftPoint_ +
					(1.0 - row) * camera_.upperLeftPoint_) +
					(1.0 - col) * (row * camera_.lowerRightPoint_ +
						(1.0 - row) * camera_.upperRightPoint_);
				Vector4d viewingRayDirection = (pointOnImagePlane - camera_.eyePoint_).normalized();
				Vector4d viewingRayPoint = camera_.eyePoint_;
				Ray viewingRay(viewingRayPoint, viewingRayDirection);

				double minDistanceClosestObj = std::numeric_limits<double>::infinity();
				// for (auto ptr : geometries_) {
				// 	Vector4d intersectionPt, normalDirection;
				// 	if (ptr -> calculateIntersectionNormal(viewingRay, intersectionPt, normalDirection) {
				// 		float distance = 
				// 	}

				// }
			}
		}
	}

	/***** CAMERA *****/
	bool hasCamera() { return hasCamera_; }
	const Camera& camera() { return camera_; }
	void camera(const Camera& camera) {
		hasCamera_ = true;
		camera_ = camera;
	}
	/***** GEOMETRIES ******/
	void addGeometry(std::unique_ptr<Geometry>&& geometry) {
		geometries_.push_back(std::move(geometry));
	}
	/***** LIGHTS ******/
	void addLight(std::unique_ptr<Light>&& light) {
		lights_.push_back(std::move(light));
	}
private:
	bool hasCamera_;
	Camera camera_;
	std::vector<std::unique_ptr<Geometry>> geometries_;
	std::vector<std::unique_ptr<Light>> lights_;
};

class RTParser {
private:
	enum LineTypeValue {
		LINE_TYPE_CAMERA,
		LINE_TYPE_SPHERE,
		LINE_TYPE_TRIANGLE,
		LINE_TYPE_POINT_LIGHT,
		LINE_TYPE_DIRECTIONAL_LIGHT,
		LINE_TYPE_AMBIENT_LIGHT,
		LINE_TYPE_MATERIAL,
		LINE_TYPE_TRANSFORM_TRANSLATE,
		LINE_TYPE_TRANSFORM_ROTATE,
		LINE_TYPE_TRANSFORM_SCALE,
		LINE_TYPE_TRANSFORM_IDENTITY
	};
	struct LineType {
		LineTypeValue type_;
		int pmin_, pmax_;
		LineType() = default;
		LineType(const LineType&) = default;
		LineType(LineTypeValue type, int pval) :
			LineType(type, pval, pval) {}
		LineType(LineTypeValue type, int pmin, int pmax):
			type_(type), pmin_(pmin), pmax_(pmax) {}
	};
	static const std::map<std::string, LineType> LINE_TYPES;
	static std::map<std::string, LineType> initializeLineTypes() {
		return {
			{"cam", {LINE_TYPE_CAMERA, 15}},
			{"sph", {LINE_TYPE_SPHERE, 4}},
			{"tri", {LINE_TYPE_TRIANGLE, 9}},
			{"ltp", {LINE_TYPE_POINT_LIGHT, 6, 7}},
			{"ltd", {LINE_TYPE_DIRECTIONAL_LIGHT, 6}},
			{"lta", {LINE_TYPE_AMBIENT_LIGHT, 3}},
			{"mat", {LINE_TYPE_MATERIAL, 13}},
			{"xft", {LINE_TYPE_TRANSFORM_TRANSLATE, 3}},
			{"xfr", {LINE_TYPE_TRANSFORM_ROTATE, 3}},
			{"xfs", {LINE_TYPE_TRANSFORM_SCALE, 3}},
			{"xfz", {LINE_TYPE_TRANSFORM_IDENTITY, 0}}
		};
	}
public:
	RTParser(GlobalScene& scene) :
		scene_(scene),
		transform_() {}
	void parseFile(std::string filename) {
		std::ifstream stream(filename);
		if (!stream)
			throw ParseException("file not found: " + filename);
		int lineno = 1;
		for (std::string line; std::getline(stream, line); lineno++) {
			std::istringstream ss(line);
			std::string stype = extractToken(ss, lineno);
			if (stype.empty()) {
				// blank line
				continue;
			}
			if (stype == "obj") {
				// obj needs special handling
				std::string filename = extractToken(ss, lineno);
				if (filename.empty())
					throw ParseException("obj requires a filename", lineno);
				// TODO
				abort();
				continue;
			}
			auto iter = LINE_TYPES.find(stype);
			if (iter == LINE_TYPES.end()) {
				std::ostringstream es;
				es << "unknown line type " << stype;
				ParseException::showWarning(es.str(), lineno);
				continue;
			}
			LineType ltype = iter->second;
			std::vector<double> params;
			for (int i = 0; i < ltype.pmax_; i++) {
				std::string token = extractToken(ss, lineno);
				if (token.empty()) {
					if (i < ltype.pmin_) {
						std::ostringstream es;
						es << stype
							<< " requires "
							<< (ltype.pmin_ == ltype.pmax_ ? "" : "at least ")
							<< ltype.pmin_
							<< " parameters";
						throw ParseException(es.str(), lineno);
					} else {
						break;
					}
				}
				try {
					params.push_back(std::stod(token));
				} catch (std::logic_error& e) {
					std::ostringstream es;
					es << "invalid number " << token;
					throw ParseException(es.str(), lineno);
				}
			}
			auto cvec = [&](int offset){
				return Color3d(&params[offset]);
			};
			auto hvec = [&](int offset) {
				return Vector3d(&params[offset]).homogeneous();
			};
			auto dvec = [&](int offset) {
				Vector3d c(&params[offset]);
				if (c.isZero())
					throw ParseException("zero direction specified", lineno);
				return Util::vec4dFrom3d(c.normalized());
			};
			Vector3d fvec = cvec(0);
			switch (ltype.type_) {
				case LINE_TYPE_TRANSFORM_IDENTITY:
					transform_.setIdentity();
					break;
				case LINE_TYPE_TRANSFORM_TRANSLATE:
					transform_.pretranslate(fvec);
					break;
				case LINE_TYPE_TRANSFORM_SCALE:
					transform_.prescale(fvec);
					break;
				case LINE_TYPE_TRANSFORM_ROTATE:
					transform_.prerotate(AngleAxis<double>(
						fvec.norm() * (2 * M_PI / 360.0),
						fvec.normalized()
					));
					break;
				case LINE_TYPE_MATERIAL:
					material_.ambientColor_ = cvec(0);
					material_.diffuseColor_ = cvec(3);
					material_.specularColor_ = cvec(6);
					material_.reflectiveColor_ = cvec(10);
					material_.specularCoefficient_ = params[9];
					break;
				case LINE_TYPE_CAMERA:
				{
					Camera camera;
					camera.transform_ = transform_;
					camera.eyePoint_ = hvec(0);
					camera.lowerLeftPoint_ = hvec(3);
					camera.lowerRightPoint_ = hvec(6);
					camera.upperLeftPoint_ = hvec(9);
					camera.upperRightPoint_ = hvec(12);
					scene_.camera(camera);
					break;
				}
				case LINE_TYPE_SPHERE:
				{
					Sphere* sphere = new Sphere();
					sphere->transform_ = transform_;
					sphere->material_ = material_;
					sphere->center_ = hvec(0);
					sphere->radius_ = params[3];
					scene_.addGeometry(std::unique_ptr<Geometry>(sphere));
					break;
				}
				case LINE_TYPE_TRIANGLE:
				{
					Mesh* mesh = new Mesh();
					mesh->transform_ = transform_;
					mesh->material_ = material_;
					mesh->addTriangle({{ hvec(0), hvec(3), hvec(6) }});
					scene_.addGeometry(std::unique_ptr<Geometry>(mesh));
					break;
				}
				case LINE_TYPE_POINT_LIGHT:
				{
					PointLight::Falloff falloff = PointLight::FALLOFF_NONE;
					switch ((int)params[6]) {
						case 0: falloff = PointLight::FALLOFF_NONE; break;
						case 1: falloff = PointLight::FALLOFF_LINEAR; break;
						case 2: falloff = PointLight::FALLOFF_QUADRATIC; break;
						default:
							ParseException::showWarning(
								"invalid falloff type", lineno);
							break;
					}
					PointLight* light = new PointLight();
					light->transform_ = transform_;
					light->point_ = hvec(0);
					light->color_ = cvec(3);
					light->falloff_ = falloff;
					scene_.addLight(std::unique_ptr<Light>(light));
					break;
				}
				case LINE_TYPE_DIRECTIONAL_LIGHT:
				{
					DirectionalLight* light = new DirectionalLight();
					light->transform_ = transform_;
					light->direction_ = dvec(0);
					light->color_ = cvec(3);
					scene_.addLight(std::unique_ptr<Light>(light));
					break;
				}
				case LINE_TYPE_AMBIENT_LIGHT:
				{
					AmbientLight* light = new AmbientLight();
					light->transform_ = transform_;
					light->color_ = cvec(0);
					scene_.addLight(std::unique_ptr<Light>(light));
					break;
				}
			}
		}
	}
	std::string extractToken(std::istream& stream, int lineno) {
		// skip beginning spaces
		int c;
		while (true) {
			c = stream.peek();
			if (c == EOF)
				return std::string();
			if (!std::isspace(c))
				break;
			stream.get();
		}
		// read string value
		std::ostringstream s;
		if (stream.peek() == '"') {
			stream.get();
			while (true) {
				c = stream.get();
				if (c == EOF)
					throw ParseException("unclosed quotes", lineno);
				if (c == '"')
					break;
				s.put(c);
			}
		} else {
			while (true) {
				c = stream.get();
				if (c == EOF || std::isspace(c))
					break;
				s.put(c);
			}
		}
		// return it
		return s.str();
	}
private:
	GlobalScene& scene_;
	Transform4d transform_;
	Material material_;
};

class PNGWriter {
private:
	typedef Matrix<uint8_t, 3, 1> Vector3uc;
	typedef Array<Vector3uc, Dynamic, Dynamic, RowMajor> RGBImage;
public:
	PNGWriter(std::string filename) : filename_(filename) {}
	void writeImage(GlobalScene::RasterImage& image) {
		RGBImage rgb = convertImage(image);
		png_image png = {0};
		png.version = PNG_IMAGE_VERSION;
		png.width = image.cols();
		png.height = image.rows();
		png.format = PNG_FORMAT_RGB;
		if (!png_image_write_to_file(&png, filename_.c_str(), false,
		                             rgb.data(), image.cols() * 3, NULL))
			throw WriteException(std::string(png.message));
	}
	RGBImage convertImage(GlobalScene::RasterImage& image) {
		RGBImage rgb(image.rows(), image.cols());
		for (int i = 0; i < image.size(); i++)
			rgb(i) = (image(i).cwiseMin(1).cwiseMax(0) * 255.0).cast<uint8_t>();
		return rgb;
	}
private:
	std::string filename_;
};

const std::map<std::string, RTParser::LineType> RTParser::LINE_TYPES(RTParser::initializeLineTypes());

static struct option programOptions[] = {
	{"help", 0, NULL, 'h'},
	{"output", 1, NULL, 'o'}
};

static void printHelp(const char *prog) {
	std::cerr << "Usage: " << prog << " [options] -o <output file> <input files>..." << std::endl;
}

int main(int argc, char *argv[]) {
	// program options
	std::string outputFilename;
	// parse options
	int opt;
	while ((opt = getopt_long(argc, argv, "ho:", programOptions, NULL)) != -1) {
		switch (opt) {
			case 'o':
				outputFilename = optarg;
				break;
			case 'h':
			case '?':
				printHelp(argv[0]);
				return 1;
		}
	}
	if (optind >= argc) {
		std::cerr << "Error: At least one input file must be specified." << std::endl;
		return 1;
	}
	if (outputFilename.empty()) {
		std::cerr << "Error: An output file must be specified." << std::endl;
		return 1;
	}
	// check that output is writable
	if (!std::ofstream(outputFilename)) {
		std::cerr << "Error: Output file is not writable." << std::endl;
		return 1;
	} else {
		std::remove(outputFilename.c_str());
	}
	// read input
	GlobalScene scene;
	while (optind < argc) {
		RTParser parser(scene);
		try {
			parser.parseFile(argv[optind++]);
		} catch (const ParseException& e) {
			std::cerr << "Error: " << e.what() << std::endl;
			return 1;
		}
	}
	if (!scene.hasCamera()) {
		std::cerr << "Error: At least one camera must be specified." << std::endl;
		return 1;
	}
	// render scene
	GlobalScene::RasterImage image(2000,3000);
	scene.renderScene(image);
	// write output
	try {
		PNGWriter writer(outputFilename);
		writer.writeImage(image);
	} catch (const WriteException& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return 1;
	}
	// done!
	return 0;
}
