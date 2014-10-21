#pragma once

class Light : public Transformable {
public:
	virtual double calculateDistanceToLight(const Vector4d& source) = 0;
	virtual Vector4d calculateDirectionToLight(const Vector4d& source) = 0;
	virtual Color3d colorForDistance(double dist) { return color_; }
	Ray calculateRayToLight(const Vector4d& source) {
		return Ray(source, calculateDirectionToLight(source));
	}
public:
	Color3d color_;
};

class PointLight : public Light {
public:
	double calculateDistanceToLight(const Vector4d& source) {
		return calculateDifferenceToLight(source).norm();
	}
	Vector4d calculateDirectionToLight(const Vector4d& source) {
		return calculateDifferenceToLight(source);
	}
	Color3d colorForDistance(double dist) {
		return pow(dist, -falloffExponent_) * color_;
	}
	void point(const Vector4d& p) { xfDirty_ = true; point_ = p; }
private:
	void updateTransformedPoint() {
		if (xfDirty_) {
			xfDirty_ = false;
			xfPoint_ = forwardTransform() * point_;
		}
	}
	Vector4d calculateDifferenceToLight(const Vector4d& source) {
		updateTransformedPoint();
		return xfPoint_ - source;
	}
private:
	Vector4d point_, xfPoint_;
	bool xfDirty_ = true;
public:
	double falloffExponent_;
};

class DirectionalLight : public Light {
public:
	double calculateDistanceToLight(const Vector4d& source) {
		return std::numeric_limits<double>::infinity();
	}
	Vector4d calculateDirectionToLight(const Vector4d& source) {
		updateTransformedDirection();
		return -xfDirection_;
	}
	void direction(const Vector4d& d) { xfDirty_ = true; direction_ = d; }
private:
	void updateTransformedDirection() {
		if (xfDirty_) {
			xfDirty_ = false;
			xfDirection_ = forwardTransform() * direction_;
		}
	}
private:
	Vector4d direction_, xfDirection_;
	bool xfDirty_ = true;
};

class AmbientLight : public Light {
public:
	double calculateDistanceToLight(const Vector4d& source) {
		return 0.0;
	}
	Vector4d calculateDirectionToLight(const Vector4d& source) {
		return Vector4d::Zero();
	}
};
