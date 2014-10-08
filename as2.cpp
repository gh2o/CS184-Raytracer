#include <getopt.h>
#include <png.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <locale>
#include <memory>
#include <map>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Geometry>

using namespace Eigen;

typedef Transform<double,3,Affine> Transform4d;
typedef Array<double,3,1> Color3d;

class Ray {
public:
	Vector4d origin_;
	Vector4d direction_;
	double t_;

	Ray(Vector4d inputOrigin, Vector4d inputDirection, double inputT) {
		origin_ = inputOrigin;
		direction_ = inputDirection;
		t_ = inputT;
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
	Transform4d objectTransform_;
};

class Camera : public Transformable {
public:
	Vector4d eyePoint_;
	Vector4d lowerLeftPoint_;
	Vector4d lowerRightPoint_;
	Vector4d upperLeftPoint_;
	Vector4d upperRightPoint_;
};

class Sphere {
private:
	Vector4d centerPt_;
	float radius_;

public:

	Sphere(Vector4d inputCenter, float inputRadius) {
		radius_ = inputRadius;
		centerPt_ = inputCenter;
	}

	Vector4d centerPt() { return centerPt_; }
	float radius() { return radius_; }

	list<Vector4d> sphereIntersection(Ray *inputRay) {
		/***** TO DO FOR LATER *****/
		//see if vector collides first; calculate projection circle center
		//onto ray and see if distance is > 0
		//Vector4d circleCtrProjection

		list<Vector4d> returnList;
		float a = inputRay->direction_.dot(inputRay->direction_);
		float b = 2 * inputRay->direction_.dot(inputRay->origin_ - centerPt_);
		float c = (inputRay->origin_ - centerPt_).dot(inputRay->origin_ - centerPt_) - radius_*radius_;
		float discriminant = pow(b,2) - 4*a*c;
		float quadraticFormResult = float(-1 * b - sqrt(b*b-4*a*c)) / float(2*a);

		Vector4d projectionPt = 
			((centerPt_.dot(inputRay->origin_)) / pow(centerPt_.norm(), 2)) * centerPt_;

		if (discriminant < 0) {
			return returnList;
		} 

		else if (discriminant == 0) {
			returnList.push_back(projectionPt);
		} else {
			Vector4d a = radius_;
			Vector4d b = projectionPt - centerPt_;
			b = b.cwiseAbs();
			Vector4d c = sqrt(pow(a,2) - b * b);
			Vector4d distRayOrigToProjPt = projectonPt - centerPt_;
			distRayOrigToProjPt = distRayOrigToProjPt.cwiseAbs();
			Vector4d distRayOrigToCenterPt = distRayOrigToProjPt - c;

			// finding first inersection point
			Vector4d intersectionPt1 = inputRay->origin_ + inputRay->direction_ * distRayOrigToCenterPt;

			// finding second intersection point
			Vector4d halfCircleSegment = projectionPt_ - intersectionPt1;
			Vector4d intersectionPt2 = intersectionPt1 + 2*halfCircleSegment;

			returnList.push_back(intersectionPt1);
			returnList.push_back(intersectionPt2);
		}
		return returnList;
	}

};

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

class WriteException : public std::runtime_error {
public:
	using runtime_error::runtime_error;
};

class GlobalScene {
public:
	typedef Array<Vector3d, Dynamic, Dynamic, RowMajor> RasterImage;
public:
	GlobalScene() :
		hasCamera_(false) {}
	void renderScene(RasterImage& output) {
		for (int i = 200; i < 600; i++)
			for (int j = 200; j < 400; j++)
				output(i,j) << 0,1,1;
	}
	bool hasCamera() { return hasCamera_; }
	const Camera& camera() { return camera_; }
	void camera(const Camera& camera) {
		hasCamera_ = true;
		camera_ = camera;
	}
private:
	bool hasCamera_;
	Camera camera_;
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
					Camera camera;
					camera.objectTransform_ = transform_;
					camera.eyePoint_ = hvec(0);
					camera.lowerLeftPoint_ = hvec(3);
					camera.lowerRightPoint_ = hvec(6);
					camera.upperLeftPoint_ = hvec(9);
					camera.upperRightPoint_ = hvec(12);
					scene_.camera(camera);
					break;
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
