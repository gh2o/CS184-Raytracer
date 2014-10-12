#include <fstream>
#include "parsers.h"
#include "util.h"

std::map<std::string, RTIParser::LineType> RTIParser::initializeLineTypes() {
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

const std::map<std::string, RTIParser::LineType> RTIParser::LINE_TYPES =
	RTIParser::initializeLineTypes();

static std::string extractToken(std::istream& stream, int lineno) {
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
	bool quoted = stream.peek() == '"';
	std::ostringstream os;
	if (quoted) {
		stream.get();
		while (true) {
			c = stream.get();
			if (c == EOF)
				throw ParseException("unclosed quotes", lineno);
			if (c == '"')
				break;
			os.put(c);
		}
	} else {
		while (true) {
			c = stream.get();
			if (c == EOF || std::isspace(c))
				break;
			os.put(c);
		}
	}
	// ignore comments
	std::string str = os.str();
	if (!quoted && !str.empty() && str[0] == '#') {
		str.clear();
		stream.ignore(std::numeric_limits<std::streamsize>::max());
	}
	// return it
	return str;
}
void RTIParser::parseFile(std::string filename) {
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
		auto cvec = [&](int offset) -> Color3d {
			return Color3d(&params[offset]);
		};
		auto tvec = [&](int offset) -> Vector3d {
			return Vector3d(&params[offset]);
		};
		auto hvec = [&](int offset) -> Vector4d {
			return Vector3d(&params[offset]).homogeneous();
		};
		auto dvec = [&](int offset) -> Vector4d {
			Vector3d c(&params[offset]);
			if (c.isZero())
				throw ParseException("zero direction specified", lineno);
			return Util::vec4dFrom3d(c.normalized());
		};
		switch (ltype.type_) {
			case LINE_TYPE_TRANSFORM_IDENTITY:
				transform_.setIdentity();
				break;
			case LINE_TYPE_TRANSFORM_TRANSLATE:
				transform_.pretranslate(tvec(0));
				break;
			case LINE_TYPE_TRANSFORM_SCALE:
				transform_.prescale(tvec(0));
				break;
			case LINE_TYPE_TRANSFORM_ROTATE:
			{
				Vector3d rvec = tvec(0);
				if (!rvec.isZero())
					transform_.prerotate(AngleAxis<double>(
						rvec.norm() * (2 * M_PI / 360.0),
						rvec.normalized()
					));
				break;
			}
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
