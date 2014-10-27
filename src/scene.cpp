#include "scene.h"
#include "options.h"
#include <cmath>

static void dummyProgressHandler(int complete, int total) {}

void Scene::renderScene(RasterImage& output, ProgressHandler phandler) {
	phandler || (phandler = dummyProgressHandler);
	for (int r = 0; r < output.rows(); r++) {
		for (int c = 0; c < output.cols(); c++) {
			phandler(r * output.cols() + c, output.size());
			double rowFrac = (r + 0.5) / output.rows();
			double colFrac = (c + 0.5) / output.cols();
			Ray viewingRay = camera_.calculateViewingRay(rowFrac, colFrac);
			output(r,c) = traceRay(viewingRay, programOptions.bounceDepth_);
		}
	}
	phandler(output.size(), output.size());
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

Color3d Scene::traceRay(Ray viewingRay, int bounceDepth, bool fromInside) {
	Geometry* targetGeometry;
	Vector4d targetIntersection, targetNormal;
	double targetDistance;
	bool rayIntersected = castRay(viewingRay, &targetDistance, &targetGeometry,
			&targetIntersection, &targetNormal, fromInside);
	if (!rayIntersected)
		return Color3d::Zero();
	if (programOptions.intersectionOnly_)
		return Color3d::Constant(1.0 / (targetDistance * targetDistance));
	// if inside, invert the targetNormal
	if (fromInside) {
		targetNormal = -targetNormal;
	}
	// normalize normal
	targetNormal.normalize();
	// calulate resulting color
	Color3d resultColor = Color3d::Zero();
	for (auto& pointer : lights_) {
		Light& light = *pointer;
		if (dynamic_cast<AmbientLight*>(&light)) {
			// ambient
			double ambientIntensity = 1.0;
			resultColor += light.color_ *
				targetGeometry->material_.ambientColor_ * ambientIntensity;
		} else {
			// check for occlusion
			Ray rayToLight = light.calculateRayToLight(targetIntersection);
			bool lightReverseNormals = targetNormal.dot(rayToLight.direction()) < 0;
			double distToLight = light.calculateDistanceToLight(targetIntersection);
			double distToOccluder;
			if (castRay(rayToLight, &distToOccluder, nullptr, nullptr, nullptr, lightReverseNormals ^ fromInside)
					&& distToOccluder <= distToLight)
				continue;
			// color
			Color3d attenuatedColor = light.colorForDistance(distToLight);
			// diffuse
			double diffuseIntensity = std::max(targetNormal.dot(rayToLight.direction()), 0.0);
			resultColor += diffuseIntensity *
				attenuatedColor * targetGeometry->material_.diffuseColor_;
			// specular
			Vector4d reflectDirection =
				2 * targetNormal.dot(rayToLight.direction()) * targetNormal - rayToLight.direction();
			double specularIntensity = pow(std::max(-viewingRay.direction().dot(reflectDirection), 0.0),
				targetGeometry->material_.specularCoefficient_ );
			resultColor += specularIntensity *
				attenuatedColor * targetGeometry->material_.specularColor_;
		}
	}
	// bounce!!!
	const Vector4d& incomingDirection = viewingRay.direction();
	const Color3d& reflectiveColor = targetGeometry->material_.reflectiveColor_;
	const Color3d& translucencyColor = targetGeometry->material_.translucencyColor_;

	if (bounceDepth > 0 && !reflectiveColor.isZero()) {
		Vector4d outgoingDirection =
			incomingDirection - 2 * targetNormal.dot(incomingDirection) * targetNormal;


		Ray outgoingRay(targetIntersection, outgoingDirection);
		resultColor += traceRay(outgoingRay, bounceDepth - 1, fromInside) * reflectiveColor;
	}

	if (bounceDepth > 0 && !translucencyColor.isZero()) {
		double indexOfRefractivity = targetGeometry->material_.indexOfRefractivity_;
		double n;
		if (fromInside) {
			n = indexOfRefractivity;
		} else {
			n = 1.0 / indexOfRefractivity;
		}
		double cosI = targetNormal.dot(incomingDirection);
		double sinT2 = n * n * (1.0 - cosI * cosI);
		if (sinT2 > 1.0) {
			// vector is invalid
		} else {
			Vector4d refractedDirection = n * incomingDirection - (n + sqrt(1.0 - sinT2)) * targetNormal;
			Ray refractedRay(targetIntersection, refractedDirection);
			resultColor += traceRay(refractedRay, bounceDepth - 1, !fromInside);
		}
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
