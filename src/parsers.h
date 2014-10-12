#pragma once
#include <map>
#include "scene.h"

class RTIParser {
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
	static std::map<std::string, LineType> initializeLineTypes();
public:
	RTIParser(Scene& scene) :
		scene_(scene),
		transform_(Transform4d::Identity()) {}
	void parseFile(std::string filename);
private:
	Scene& scene_;
	Transform4d transform_;
	Material material_;
};
