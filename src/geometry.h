#pragma once
#include <array>
#include <vector>
#include "rtbase.h"

class Geometry : public Transformable {
public:
	bool calculateIntersectionNormal(Ray inputRay, Vector4d& intersectionPt, Vector4d& normalDirection,
			bool reverseNormals);
	virtual bool calculateIntNormInObjSpace(Ray inputRay, Vector4d& intersectionPt, Vector4d& normalDirection,
			bool reverseNormals) = 0;
public:
	Material material_;
};

class Sphere : public Geometry {
public:
	bool calculateIntNormInObjSpace(Ray inputRay, Vector4d& intersectionPt, Vector4d& normalDirection,
			bool reverseNormals);
public:
	Vector4d center_;
	float radius_;
};

class Mesh : public Geometry {
public:
	void addTriangle(const std::array<Vector4d,3>& points);
	bool calculateIntNormInObjSpace(Ray inputRay, Vector4d& intersectionPt, Vector4d& normalDirection,
			bool reverseNormals);
public:
	struct Face { std::array<Vector4d,3> points_, normals_; };
	std::vector<Face> faces_;
};
