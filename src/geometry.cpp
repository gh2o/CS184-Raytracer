#include "geometry.h"
#include "util.h"

bool Geometry::calculateIntersectionNormal(Ray inputRay, Vector4d& intersectionPt, Vector4d& normalDirection,
		bool reverseNormals) {
	Matrix4d fwdTransform = transform_.matrix();
	Matrix4d invTransform = transform_.inverse().matrix();
	Ray transformedRay(invTransform * inputRay.origin(), invTransform * inputRay.direction());
	Vector4d transformedIntersectionPt, transformedNormalDirection;
	bool hasIntersection = calculateIntNormInObjSpace(transformedRay, transformedIntersectionPt, 
		transformedNormalDirection, reverseNormals);
	if (!hasIntersection) {
		return false;
	}
	intersectionPt = fwdTransform * transformedIntersectionPt;
	normalDirection = invTransform.transpose() * transformedNormalDirection;
	normalDirection.w() = 0; // translation part of matrix introduces residue
	if (fwdTransform.determinant() < 0) {
		normalDirection = -normalDirection;
	}
	return true;
}

bool Sphere::calculateIntNormInObjSpace(Ray inputRay, Vector4d& intersectionPt, Vector4d& normalDirection,
		bool reverseNormals) {
	// always assume ray coming from outside of sphere
	Vector4d originCenterDiff = inputRay.origin() - center_;
	double a = inputRay.direction().dot(inputRay.direction());
	double b = 2 * inputRay.direction().dot(originCenterDiff);
	double c = originCenterDiff.dot(originCenterDiff) - radius_*radius_;
	double discriminant = b*b - 4*a*c;
	double result;
	if (discriminant < 0)
		return false;
	if (reverseNormals)
		result = (-b + sqrt(discriminant)) / (2 * a);
	else
		result = (-b - sqrt(discriminant)) / (2 * a);
	if (result < 0)
		return false;
	intersectionPt = inputRay.origin() + result * inputRay.direction();
	normalDirection = intersectionPt - center_;
	return true;
}

bool Mesh::calculateIntNormInObjSpace(Ray inputRay, Vector4d& intersectionPt, Vector4d& normalDirection,
		bool reverseNormals) {
	bool closestExists = false;
	double closestDist = std::numeric_limits<double>::infinity();
	for (Face& face : faces_) {
		// base vectors for matrix
		Vector4d va = face.points_[1] - face.points_[0];
		Vector4d vb = face.points_[2] - face.points_[0];
		Vector4d dir = inputRay.direction();
		Matrix4d m;
		m << va, vb, -dir, Vector4d(0,0,0,1);
		if (m.determinant() == 0)
			continue;
		// solve matrix for va and vb coeffs, direction scale t
		Vector4d sol = m.inverse() * (inputRay.origin() - face.points_[0]);
		double a = sol(0);
		double b = sol(1);
		double t = sol(2);
		// validate constraints
		if (a < 0 || b < 0 || a + b > 1)
			continue;
		if (t < 0)
			continue;
		// check if closest
		double dist = t * dir.norm();
		if (dist >= closestDist)
			continue;
		// check normal if one sided
		Vector4d testNormal =
			(1.0 - a - b) * face.normals_[0] +
			a * face.normals_[1] +
			b * face.normals_[2];
		bool hitsFront = testNormal.dot(inputRay.direction()) < 0;
		if (!hitsFront ^ reverseNormals)
			continue; // normal should be opposite of ray direction
		// update closest
		closestExists = true;
		closestDist = dist;
		intersectionPt = face.points_[0] + a * va + b * vb;
		normalDirection = testNormal;
	}
	return closestExists;
}

void Mesh::addTriangle(const std::array<Vector4d,3>& points) {
	const Vector4d& v0 = points[0];
	const Vector4d& v1 = points[1];
	const Vector4d& v2 = points[2];
	Vector4d n = Util::cross(v1 - v0, v2 - v0).normalized();
	Vector4d ep = std::numeric_limits<double>::epsilon()
		* (v0 + v1 + v2).norm() / 3 * n;
	for (int s = -1; s <= 1; s += 2) {
		Face face;
		face.points_ = points;
		face.normals_ = {{ s*n, s*n, s*n }};
		for (Vector4d& point : face.points_)
			point += s*ep;
		faces_.push_back(face);
	}
}
