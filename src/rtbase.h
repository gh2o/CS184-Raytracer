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
	Color3d translucencyColor_;
	double indexOfRefractivity_;
};

class Transformable {
public:
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
	Transform4d forwardTransform_ = Transform4d::Identity();
	Transform4d inverseTransform_ = Transform4d::Identity();
	double transformDeterminant_;
};

class Camera : public Transformable {
public:
	void eyePoint(const Vector4d& p)        { xfDirty_ = true; eyePoint_ = p; }
	void lowerLeftPoint(const Vector4d& p)  { xfDirty_ = true; lowerLeftPoint_ = p; }
	void lowerRightPoint(const Vector4d& p) { xfDirty_ = true; lowerRightPoint_ = p; }
	void upperLeftPoint(const Vector4d& p)  { xfDirty_ = true; upperLeftPoint_ = p; }
	void upperRightPoint(const Vector4d& p) { xfDirty_ = true; upperRightPoint_ = p; }
public:
	Ray calculateViewingRay(double rowFrac, double colFrac) {
		updateTransformedVectors();
		Vector4d imagePlanePoint =
			colFrac * (
				rowFrac *         xfLowerRightPoint_ +
				(1.0 - rowFrac) * xfUpperRightPoint_) +
			(1.0 - colFrac) * (
				rowFrac         * xfLowerLeftPoint_ +
				(1.0 - rowFrac) * xfUpperLeftPoint_);
		return Ray(xfEyePoint_, imagePlanePoint - xfEyePoint_);
	}
private:
	void updateTransformedVectors() {
		if (xfDirty_) {
			xfDirty_ = false;
			xfEyePoint_        = forwardTransform() * eyePoint_;
			xfLowerLeftPoint_  = forwardTransform() * lowerLeftPoint_;
			xfLowerRightPoint_ = forwardTransform() * lowerRightPoint_;
			xfUpperLeftPoint_  = forwardTransform() * upperLeftPoint_;
			xfUpperRightPoint_ = forwardTransform() * upperRightPoint_;
		}
	}
private:
	Vector4d eyePoint_,        xfEyePoint_;
	Vector4d lowerLeftPoint_,  xfLowerLeftPoint_;
	Vector4d lowerRightPoint_, xfLowerRightPoint_;
	Vector4d upperLeftPoint_,  xfUpperLeftPoint_;
	Vector4d upperRightPoint_, xfUpperRightPoint_;
	bool xfDirty_ = true;
};
