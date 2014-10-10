#include "scene.h"
#include "options.h"

void Scene::renderScene(RasterImage& output) {
	for (int r = 0; r < output.rows(); r++) {
		for (int c = 0; c < output.cols(); c++) {
			double row = (r + 0.5) / output.rows();
			double col = (c + 0.5) / output.cols();
			Vector4d pointOnImagePlane =
				col * (
					row *         camera_.lowerRightPoint_ +
					(1.0 - row) * camera_.upperRightPoint_) +
				(1.0 - col) * (
					row         * camera_.lowerLeftPoint_ +
					(1.0 - row) * camera_.upperLeftPoint_);
			Ray viewingRay(camera_.eyePoint_, pointOnImagePlane - camera_.eyePoint_);
			output(r,c) = traceRay(viewingRay, programOptions.bounceDepth_);
		}
	}
	if (programOptions.intersectionOnly_) {
		double maxBrightness = std::numeric_limits<double>::min();
		for (int r = 0; r < output.rows(); r++)
			for (int c = 0; c < output.cols(); c++)
				maxBrightness = std::max(maxBrightness, output(r,c).maxCoeff());
		for (int r = 0; r < output.rows(); r++)
			for (int c = 0; c < output.cols(); c++)
				output(r,c) /= maxBrightness;
	}
}

Color3d Scene::traceRay(Ray viewingRay, int bounceDepth) {
	Geometry* targetGeometry;
	Vector4d targetIntersection, targetNormal;
	double targetDistance;
	bool rayIntersected = castRay(viewingRay, &targetDistance, &targetGeometry,
			&targetIntersection, &targetNormal);
	if (!rayIntersected)
		return Color3d::Zero();
	if (programOptions.intersectionOnly_)
		return Color3d::Constant(1.0 / (targetDistance * targetDistance));
	// normalize normal
	targetNormal.normalize();
	// calulate resulting color
	Color3d resultColor = Color3d::Zero();
	for (auto& pointer : lights_) {
		Light& light = *pointer;
		if (dynamic_cast<AmbientLight*>(&light)) {
			// ambient
			double ambientIntensity = 1.0;
			resultColor += ambientIntensity * light.color_ *
				targetGeometry->material_.ambientColor_;
		} else {
			// check for occlusion
			Ray rayToLight = light.calculateRayToLight(targetIntersection);
			double distToLight = light.calculateDistanceToLight(targetIntersection);
			double distToOccluder;
			if (castRay(rayToLight, &distToOccluder, nullptr, nullptr, nullptr)
					&& distToOccluder <= distToLight)
				continue;
			// diffuse
			double diffuseIntensity = std::max(targetNormal.dot(rayToLight.direction()), 0.0);
			resultColor += diffuseIntensity *
				light.color_ * targetGeometry->material_.diffuseColor_;
			// specular
			Vector4d reflectDirection =
				2 * targetNormal.dot(rayToLight.direction()) * targetNormal - rayToLight.direction();
			double specularIntensity = pow(std::max(-viewingRay.direction().dot(reflectDirection), 0.0),
				targetGeometry->material_.specularCoefficient_ );
			resultColor += specularIntensity *
				light.color_ * targetGeometry->material_.specularColor_;
		}
	}
	// bounce!!!
	const Color3d& reflectiveColor = targetGeometry->material_.reflectiveColor_;
	if (bounceDepth > 0 && !reflectiveColor.isZero()) {
		const Vector4d& incomingDirection = viewingRay.direction();
		Vector4d outgoingDirection =
			incomingDirection - 2 * targetNormal.dot(incomingDirection) * targetNormal;
		Ray outgoingRay(targetIntersection, outgoingDirection);
		resultColor += traceRay(outgoingRay, bounceDepth - 1) * reflectiveColor;
	}
	// done!!!
	return resultColor;
}

bool Scene::castRay(Ray castedRay, double* targetDistance, Geometry** targetGeometry,
		Vector4d* targetIntersection, Vector4d* targetNormal, bool reverseNormals) {
	double tmpDistance;
	targetDistance || (targetDistance = &tmpDistance);
	bool rayIntersected = false;
	for (auto& pointer : geometries_) {
		Geometry& testGeometry = *pointer;
		Vector4d testIntersection, testNormal;
		if (!testGeometry.calculateIntersectionNormal(castedRay, testIntersection, testNormal, reverseNormals))
			continue;
		// check if closest
		double testDistance = (testIntersection - castedRay.origin()).norm();
		if (rayIntersected && testDistance >= *targetDistance)
			continue;
		// update closest
		rayIntersected = true;
		*targetDistance = testDistance;
		if (targetGeometry)
			*targetGeometry = &testGeometry;
		if (targetIntersection)
			*targetIntersection = testIntersection;
		if (targetNormal)
			*targetNormal = testNormal;
	}
	return rayIntersected;
}
