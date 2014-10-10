#pragma once

class Light : public Transformable {
public:
	virtual double calculateDistanceToLight(const Vector4d& source) = 0;
	virtual Vector4d calculateDirectionToLight(const Vector4d& source) = 0;
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
public:
	enum Falloff {
		FALLOFF_NONE,
		FALLOFF_LINEAR,
		FALLOFF_QUADRATIC,
	};
	Vector4d point_;
	Falloff falloff_;
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
