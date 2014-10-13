#pragma once
#include "common.h"
#include "exceptions.h"

class Ray {
public:
	Ray(Vector4d ori, Vector4d dir) {
		origin(ori);
		direction(dir);
	}
	const Vector4d& origin() { return origin_; }
	const Vector4d& direction() { return direction_; }
	void origin(const Vector4d& ori) {
		if (ori.w() == 0)
			throw MathException("ray origin is a direction vector");
		origin_ = ori;
	}
	void direction(const Vector4d& dir) {
		if (dir.isZero())
			throw MathException("ray has no direction");
		if (dir.w() != 0)
			throw MathException("ray direction is a point vector");
		direction_ = dir.normalized();
	}
private:
	Vector4d origin_;
	Vector4d direction_;
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
	Transformable() :
		forwardTransform_(Transform4d::Identity()),
		inverseTransform_(Transform4d::Identity()) {}
	const Transform4d& forwardTransform() { return forwardTransform_; }
	const Transform4d& inverseTransform() { return inverseTransform_; }
	double transformDeterminant() { return transformDeterminant_; }
	void forwardTransform(const Transform4d& xf) {
		forwardTransform_ = xf;
		inverseTransform_ = xf.inverse();
		updateDeterminant();
	}
	void inverseTransform(const Transform4d& xf) {
		forwardTransform_ = xf.inverse();
		inverseTransform_ = xf;
		updateDeterminant();
	}
private:
	void updateDeterminant() {
		transformDeterminant_ = forwardTransform_.matrix().determinant();
	}
private:
	Transform4d forwardTransform_;
	Transform4d inverseTransform_;
	double transformDeterminant_;
};

class Camera : public Transformable {
public:
	Vector4d eyePoint_;
	Vector4d lowerLeftPoint_;
	Vector4d lowerRightPoint_;
	Vector4d upperLeftPoint_;
	Vector4d upperRightPoint_;
};
