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

static std::vector<std::string> extractTokens(std::istream& stream, int lineno) {
	std::vector<std::string> tokens;
	while (true) {
		std::string token = extractToken(stream, lineno);
		if (token.empty())
			break;
		tokens.push_back(token);
	}
	return tokens;
}

static std::vector<double> extractDoubles(std::istream& stream, int lineno) {
	std::vector<std::string> tokens = extractTokens(stream, lineno);
	std::vector<double> doubles(tokens.size());
	std::transform(tokens.begin(), tokens.end(), doubles.begin(),
			[&](const std::string& s){
		try {
			return std::stod(s);
		} catch (std::logic_error& e) {
			std::ostringstream es;
			es << "invalid number " << s;
			throw ParseException(es.str(), lineno);
		}
	});
	return doubles;
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
			std::string objfname = extractToken(ss, lineno);
			if (objfname.empty())
				throw ParseException("obj requires a filename", lineno);
			// convert relative path to .rti relative path
			if (objfname[0] != '/')
				objfname = Util::dirname(filename) + "/" + objfname;
			// parse the obj file
			Mesh* mesh = new Mesh();
			mesh->forwardTransform(transform_);
			mesh->material_ = material_;
			OBJParser objp(*mesh);
			objp.parseFile(objfname);
			scene_.addGeometry(std::unique_ptr<Geometry>(mesh));
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
		std::vector<double> params = extractDoubles(ss, lineno);
		if (params.size() < ltype.pmin_) {
			std::ostringstream es;
			es << stype
				<< " requires "
				<< (ltype.pmin_ == ltype.pmax_ ? "" : "at least ")
				<< ltype.pmin_
				<< " parameters";
			throw ParseException(es.str(), lineno);
		} else if (params.size() > ltype.pmax_) {
			ParseException::showWarning("extra parameters found", lineno);
		}
		while (params.size() < ltype.pmax_) {
			params.push_back(0.0);
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
				camera.forwardTransform(transform_);
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
				sphere->forwardTransform(transform_);
				sphere->material_ = material_;
				sphere->center_ = hvec(0);
				sphere->radius_ = params[3];
				scene_.addGeometry(std::unique_ptr<Geometry>(sphere));
				break;
			}
			case LINE_TYPE_TRIANGLE:
			{
				Mesh* mesh = new Mesh();
				mesh->forwardTransform(transform_);
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
				light->forwardTransform(transform_);
				light->point_ = hvec(0);
				light->color_ = cvec(3);
				light->falloff_ = falloff;
				scene_.addLight(std::unique_ptr<Light>(light));
				break;
			}
			case LINE_TYPE_DIRECTIONAL_LIGHT:
			{
				DirectionalLight* light = new DirectionalLight();
				light->forwardTransform(transform_);
				light->direction_ = dvec(0);
				light->color_ = cvec(3);
				scene_.addLight(std::unique_ptr<Light>(light));
				break;
			}
			case LINE_TYPE_AMBIENT_LIGHT:
			{
				AmbientLight* light = new AmbientLight();
				light->forwardTransform(transform_);
				light->color_ = cvec(0);
				scene_.addLight(std::unique_ptr<Light>(light));
				break;
			}
		}
	}
}

void OBJParser::parseFile(std::string filename) {
	std::ifstream stream(filename);
	if (!stream)
		throw ParseException("file not found: " + filename);
	std::vector<Vector4d> vertices;
	std::vector<Vector4d> normals;
	vertices.emplace_back(); // 1-indexed
	normals.emplace_back();  // 1-indexed
	int lineno = 1;
	for (std::string line; std::getline(stream, line); lineno++) {
		std::istringstream ss(line);
		std::string ltype = extractToken(ss, lineno);
		if (ltype.empty()) {
			// blank line
			continue;
		}
		if (ltype == "f") {
			// face
			struct OBJPoint {
				union {
					struct { int vertexIndex_, textureIndex_, normalIndex_; };
					int allIndices_[3];
				};
				Vector4d *vertexPtr_, *normalPtr_;
			};
			auto tokens = extractTokens(ss, lineno);
			if (tokens.size() < 3)
				throw ParseException("f requires at least 3 vertices", lineno);
			std::vector<OBJPoint> points(tokens.size());
			std::transform(tokens.begin(), tokens.end(), points.begin(),
					[&](const std::string& s) {
				OBJPoint p;
				const int kUnknownIndex = 0;
				int indexCount = 0;
				int prevDelimiterEnd = 0;
				while (indexCount < 3 && prevDelimiterEnd < s.length()) {
					int nextDelimiterStart = s.find("/", prevDelimiterEnd);
					if (nextDelimiterStart == std::string::npos)
						nextDelimiterStart = s.length();
					std::string currentPart = s.substr(
								prevDelimiterEnd, nextDelimiterStart);
					prevDelimiterEnd = nextDelimiterStart + 1;
					int indexValue = kUnknownIndex;
					if (!currentPart.empty()) {
						try {
							indexValue = std::stoi(currentPart);
						} catch (std::logic_error& e) {
							std::ostringstream os;
							os << "invalid integer " << currentPart;
							throw ParseException(os.str(), lineno);
						}
						if (indexValue <= 0)
							throw ParseException("index must be positive", lineno);
					}
					p.allIndices_[indexCount++] = indexValue;
				}
				while (indexCount < 3)
					p.allIndices_[indexCount++] = kUnknownIndex;
				// assign vertex
				if (p.vertexIndex_ == kUnknownIndex)
					throw ParseException("vertex index is required", lineno);
				else if (p.vertexIndex_ < vertices.size())
					p.vertexPtr_ = &vertices[p.vertexIndex_];
				else
					throw ParseException("vertex index out of range", lineno);
				// assign normal
				if (p.normalIndex_ == kUnknownIndex)
					p.normalPtr_ = nullptr;
				else if (p.normalIndex_ < normals.size())
					p.normalPtr_ = &normals[p.normalIndex_];
				else
					throw ParseException("normal index out of range", lineno);
				// return the point
				return p;
			});
			// add triangles
			OBJPoint& basePoint = points[0];
			for (auto it = points.begin() + 1; it < points.end() - 1; it++) {
				OBJPoint& firstPoint = *it;
				OBJPoint& secondPoint = *(it + 1);
				Vector4d firstVec = *firstPoint.vertexPtr_ - *basePoint.vertexPtr_;
				Vector4d secondVec = *secondPoint.vertexPtr_ - *basePoint.vertexPtr_;
				Vector4d calculatedNormal = Util::cross(firstVec, secondVec);
				if (calculatedNormal.isZero())
					throw ParseException("degenerate face", lineno);
				calculatedNormal.normalize();
				// add face to mesh
				Mesh::Face triFace;
				OBJPoint *triPoints[3] = { &basePoint, &firstPoint, &secondPoint };
				for (int i = 0; i < 3; i++) {
					triFace.points_[i] = *triPoints[i]->vertexPtr_;
					triFace.normals_[i] = triPoints[i]->normalPtr_ ?
						*triPoints[i]->normalPtr_ : calculatedNormal;
				}
				mesh_.faces_.push_back(triFace);
			}
		} else if (ltype == "v") {
			// vertex
			auto params = extractDoubles(ss, lineno);
			if (params.size() != 3 && params.size() != 4)
				throw ParseException("v requires 3 or 4 parameters", lineno);
			params.push_back(1.0); // w
			Vector4d vtx(&params[0]);
			if (vtx.w() == 0)
				throw ParseException("v must be a point vector", lineno);
			vertices.push_back(vtx);
		} else if (ltype == "vn") {
			// vertex normal
			auto params = extractDoubles(ss, lineno);
			if (params.size() != 3)
				throw ParseException("vn requires 3 parameters", lineno);
			params.push_back(0.0);
			Vector4d nrm(&params[0]);
			normals.push_back(nrm);
		} else {
			ParseException::showWarning("unknown obj line type " + ltype,
					lineno);
		}
	}
}
