#pragma once
#include "rtbase.h"

class Geometry : public Transformable {
public:
	bool calculateIntersectionNormal(Ray inputRay, Vector4d& intersectionPt, Vector4d& normalDirection);
	virtual bool calculateIntNormInObjSpace(Ray inputRay, Vector4d& intersectionPt, Vector4d& normalDirection) = 0;
public:
	Material material_;
};

class Sphere : public Geometry {
public:
	bool calculateIntNormInObjSpace(Ray inputRay, Vector4d& intersectionPt, Vector4d& normalDirection);
public:
	Vector4d center_;
	float radius_;
};

class Mesh : public Geometry {
public:
	void addTriangle(const std::array<Vector4d,3>& points);
	bool calculateIntNormInObjSpace(Ray inputRay, Vector4d& intersectionPt, Vector4d& normalDirection);
public:
	struct Face { std::array<Vector4d,3> points_, normals_; };
	std::vector<Face> faces_;
};
