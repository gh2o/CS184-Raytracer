#pragma once
#include <memory>
#include <vector>
#include "common.h"
#include "rtbase.h"
#include "geometry.h"
#include "lights.h"

class Scene {
public:
	typedef Array<Color3d, Dynamic, Dynamic, RowMajor> RasterImage;
public:
	Scene() :
		hasCamera_(false) {}
	void renderScene(RasterImage& output);
	Color3d traceRay(Ray viewingRay, int bounceDepth);
	bool castRay(Ray castedRay, double* targetDistance, Geometry** targetGeometry,
			Vector4d* targetIntersection, Vector4d* targetNormal, bool reverseNormals = false);

	/***** CAMERA *****/
	bool hasCamera() { return hasCamera_; }
	const Camera& camera() { return camera_; }
	void camera(const Camera& camera) {
		hasCamera_ = true;
		camera_ = camera;
	}
	/***** GEOMETRIES ******/
	void addGeometry(std::unique_ptr<Geometry>&& geometry) {
		geometries_.push_back(std::move(geometry));
	}
	/***** LIGHTS ******/
	void addLight(std::unique_ptr<Light>&& light) {
		lights_.push_back(std::move(light));
	}
private:
	bool hasCamera_;
	Camera camera_;
	std::vector<std::unique_ptr<Geometry>> geometries_;
	std::vector<std::unique_ptr<Light>> lights_;
};
