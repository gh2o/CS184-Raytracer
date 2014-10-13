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
		return (point_ - source).norm();
	}
	Vector4d calculateDirectionToLight(const Vector4d& source) {
		return point_ - source;
	}
	Color3d colorForDistance(double dist) {
		return pow(dist, -falloffExponent_) * color_;
	}
public:
	Vector4d point_;
	double falloffExponent_;
};

class DirectionalLight : public Light {
public:
	double calculateDistanceToLight(const Vector4d& source) {
		return std::numeric_limits<double>::infinity();
	}
	Vector4d calculateDirectionToLight(const Vector4d& source) {
		return -direction_;
	}
public:
	Vector4d direction_;
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
