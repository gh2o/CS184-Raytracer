#include <unordered_set>
#include "geometry.h"
#include "util.h"

static bool hitsBoundingBox(Ray ray, const Vector4d& bbmin, const Vector4d& bbmax) {
	for (int axis = 0; axis < 3; axis++) {
		for (int bn = 0; bn < 2; bn++) {
			double axisMag = ray.direction()[axis];
			if (axisMag == 0.0)
				continue;
			double t = ((bn ? bbmax : bbmin)[axis] - ray.origin()[axis]) / axisMag;
			if (t < 0)
				continue;
			Vector4d intersection = ray.origin() + t * ray.direction();
			for (int axis2 = 0; axis2 < 3; axis2++) {
				if (axis2 == axis)
					continue;
				if (intersection[axis2] < bbmin[axis2])
					goto missedSide;
				if (intersection[axis2] > bbmax[axis2])
					goto missedSide;
			}
			return true;
			missedSide:
			continue;
		}
	}
	return false;
}

bool Geometry::calculateIntersectionNormal(Ray inputRay, Vector4d& intersectionPt, Vector4d& normalDirection,
		bool reverseNormals) {
	Ray transformedRay(inverseTransform() * inputRay.origin(), inverseTransform() * inputRay.direction());
	Vector4d transformedIntersectionPt, transformedNormalDirection;
	bool hasIntersection = calculateIntNormInObjSpace(transformedRay, transformedIntersectionPt, 
		transformedNormalDirection, reverseNormals);
	if (!hasIntersection)
		return false;
	intersectionPt = forwardTransform() * transformedIntersectionPt;
	normalDirection = inverseTransform().matrix().transpose() * transformedNormalDirection;
	normalDirection.w() = 0; // translation part of matrix introduces residue
	if (transformDeterminant() < 0)
		normalDirection = -normalDirection;
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
	// do bounding box
	if (boundingBoxMin_ != boundingBoxMax_ && faces_.size() > 1)
		if (!hitsBoundingBox(inputRay, boundingBoxMin_, boundingBoxMax_))
			return false;
	// calculate closest face
	bool closestExists = false;
	double closestDist = std::numeric_limits<double>::infinity();
	for (Face& face : faces_) {
		// base vectors for matrix
		Vector3d va = Util::vec3dFrom4d(face.points_[1] - face.points_[0]);
		Vector3d vb = Util::vec3dFrom4d(face.points_[2] - face.points_[0]);
		Vector3d dir = Util::vec3dFrom4d(inputRay.direction());
		Vector3d rhs = Util::vec3dFrom4d(inputRay.origin() - face.points_[0]);
		Matrix3d mlower, mupper;
		mlower << va, vb, -dir;
		double dlower = mlower.determinant();
		if (dlower == 0)
			continue;
		// solve for a and ensure it's in the range
		mupper = mlower;
		mupper.col(0) = rhs;
		double a = mupper.determinant() / dlower;
		if (a < 0 || a > 1)
			continue;
		// solve for b and ensure it's in the range
		mupper = mlower;
		mupper.col(1) = rhs;
		double b = mupper.determinant() / dlower;
		if (b < 0 || a + b > 1)
			continue;
		// solve for t and ensure it's positive
		mupper = mlower;
		mupper.col(2) = rhs;
		double t = mupper.determinant() / dlower; // cramer's rule
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
		intersectionPt = face.points_[0] + Util::vec4dFrom3d(a * va + b * vb);
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

void Mesh::updateBoundingBox() {
	if (faces_.empty()) {
		boundingBoxMin_ = Vector4d::Zero();
		boundingBoxMax_ = Vector4d::Zero();
		return;
	}
	const double inf = std::numeric_limits<double>::infinity();
	boundingBoxMin_.fill(inf);
	boundingBoxMax_.fill(-inf);
	for (const Face& face : faces_) {
		for (const Vector4d& vtx : face.points_) {
			boundingBoxMin_ = boundingBoxMin_.cwiseMin(vtx);
			boundingBoxMax_ = boundingBoxMax_.cwiseMax(vtx);
		}
	}
	if (boundingBoxMin_.w() != 1.0 || boundingBoxMax_.w() != 1.0)
		throw MathException("non-unity-homogeneous bounding box");
}
