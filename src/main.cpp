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
#include "common.h"
#include "exceptions.h"
#include "options.h"
#include "util.h"
#include "rtbase.h"
#include "lights.h"
#include "geometry.h"

class GlobalScene {
public:
	typedef Array<Color3d, Dynamic, Dynamic, RowMajor> RasterImage;
public:
	GlobalScene() :
		hasCamera_(false) {}
	void renderScene(RasterImage& output) {
		for (int r = 0; r < output.rows(); r++) {
			for (int c = 0; c < output.cols(); c++) {
				double row = (r + 0.5) / output.rows();
				double col = (c + 0.5) / output.cols();
				Vector4d pointOnImagePlane =
					col * (
						row *         camera_.lowerRightPoint_ +
						(1.0 - row) * camera_.upperRightPoint_) +
					(1.0 - col) * (
						row         * camera_.lowerLeftPoint_ +
						(1.0 - row) * camera_.upperLeftPoint_);
				Ray viewingRay(camera_.eyePoint_, pointOnImagePlane - camera_.eyePoint_);
				output(r,c) = traceRay(viewingRay, programOptions.bounceDepth_);
			}
		}
		if (programOptions.intersectionOnly_) {
			double maxBrightness = std::numeric_limits<double>::min();
			for (int r = 0; r < output.rows(); r++)
				for (int c = 0; c < output.cols(); c++)
					maxBrightness = std::max(maxBrightness, output(r,c).maxCoeff());
			for (int r = 0; r < output.rows(); r++)
				for (int c = 0; c < output.cols(); c++)
					output(r,c) /= maxBrightness;
		}
	}

	Color3d traceRay(Ray viewingRay, int bounceDepth) {
		Geometry* targetGeometry;
		Vector4d targetIntersection, targetNormal;
		double targetDistance;
		bool rayIntersected = castRay(viewingRay, &targetDistance, &targetGeometry,
				&targetIntersection, &targetNormal);
		if (!rayIntersected)
			return Color3d::Zero();
		if (programOptions.intersectionOnly_)
			return Color3d::Constant(1.0 / (targetDistance * targetDistance));
		// normalize normal
		targetNormal.normalize();
		// calulate resulting color
		Color3d resultColor = Color3d::Zero();
		for (auto& pointer : lights_) {
			Light& light = *pointer;
			if (dynamic_cast<AmbientLight*>(&light)) {
				// ambient
				double ambientIntensity = 1.0;
				resultColor += ambientIntensity * light.color_ *
					targetGeometry->material_.ambientColor_;
			} else {
				// check for occlusion
				Ray rayToLight = light.calculateRayToLight(targetIntersection);
				double distToLight = light.calculateDistanceToLight(targetIntersection);
				double distToOccluder;
				if (castRay(rayToLight, &distToOccluder, nullptr, nullptr, nullptr)
						&& distToOccluder <= distToLight)
					continue;
				// diffuse
				double diffuseIntensity = std::max(targetNormal.dot(rayToLight.direction()), 0.0);
				resultColor += diffuseIntensity *
					light.color_ * targetGeometry->material_.diffuseColor_;
				// specular
				Vector4d reflectDirection =
					2 * targetNormal.dot(rayToLight.direction()) * targetNormal - rayToLight.direction();
				double specularIntensity = pow(std::max(-viewingRay.direction().dot(reflectDirection), 0.0),
					targetGeometry->material_.specularCoefficient_ );
				resultColor += specularIntensity *
					light.color_ * targetGeometry->material_.specularColor_;
			}
		}
		// bounce!!!
		const Color3d& reflectiveColor = targetGeometry->material_.reflectiveColor_;
		if (bounceDepth > 0 && !reflectiveColor.isZero()) {
			const Vector4d& incomingDirection = viewingRay.direction();
			Vector4d outgoingDirection =
				incomingDirection - 2 * targetNormal.dot(incomingDirection) * targetNormal;
			Ray outgoingRay(targetIntersection, outgoingDirection);
			resultColor += traceRay(outgoingRay, bounceDepth - 1) * reflectiveColor;
		}
		// done!!!
		return resultColor;
	}

	bool castRay(Ray castedRay, double* targetDistance, Geometry** targetGeometry,
			Vector4d* targetIntersection, Vector4d* targetNormal) {
		double tmpDistance;
		targetDistance || (targetDistance = &tmpDistance);
		bool rayIntersected = false;
		for (auto& pointer : geometries_) {
			Geometry& testGeometry = *pointer;
			Vector4d testIntersection, testNormal;
			if (!testGeometry.calculateIntersectionNormal(castedRay, testIntersection, testNormal))
				continue;
			// check if closest
			double testDistance = (testIntersection - castedRay.origin()).norm();
			if (rayIntersected && testDistance >= *targetDistance)
				continue;
			// update closest
			rayIntersected = true;
			*targetDistance = testDistance;
			if (targetGeometry)
				*targetGeometry = &testGeometry;
			if (targetIntersection)
				*targetIntersection = testIntersection;
			if (targetNormal)
				*targetNormal = testNormal;
		}
		return rayIntersected;
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
		transform_(Transform4d::Identity()) {}
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
						params.push_back(0.0);
						continue;
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
			if (!extractToken(ss, lineno).empty()) {
				ParseException::showWarning("extra parameters found", lineno);
			}
			auto cvec = [&](int offset){
				return Color3d(&params[offset]);
			};
			auto hvec = [&](int offset) -> Vector4d {
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
					if (!fvec.isZero())
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

const std::map<std::string, RTParser::LineType> RTParser::LINE_TYPES(RTParser::initializeLineTypes());

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

int main(int argc, char *argv[]) {
	// program options
	if (!programOptions.parseCommandLine(argc, argv))
		return 1;
	// check that output is writable
	const std::string& outputFilename = programOptions.outputFilename_;
	if (!std::ofstream(outputFilename)) {
		std::cerr << "Error: Output file is not writable." << std::endl;
		return 1;
	} else {
		std::remove(outputFilename.c_str());
	}
	// read input
	GlobalScene scene;
	for (std::string inputFilename : programOptions.inputFilenames_) {
		RTParser parser(scene);
		try {
			parser.parseFile(inputFilename);
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
	GlobalScene::RasterImage image(programOptions.renderHeight_,
	                               programOptions.renderWidth_);
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
