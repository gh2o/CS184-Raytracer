#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <list>


using namespace Eigen;
using namespace std;

typedef Matrix<double, 3, 1> Color;

struct {
	Color ambientColor;
	Color specularColor;
	Color reflectiveColor;
} Material;

class Ray {
public:
	Vector4d origin_;
	Vector4d direction_;
	double t_;

	Ray(Vector4d inputOrigin, Vector4d inputDirection, double inputT) {
		origin_ = inputOrigin;
		direction_ = inputDirection;
		t_ = inputT;
	}
};

class Sphere {
private:
	Vector4d centerPt_;
	float radius_;

public:

	Sphere(Vector4d inputCenter, float inputRadius) {
		radius_ = inputRadius;
		centerPt_ = inputCenter;
	}

	Vector4d centerPt() { return centerPt_; }
	float radius() { return radius_; }

	list<Vector4d> sphereIntersection(Ray *inputRay) {
		/***** TO DO FOR LATER *****/
		//see if vector collides first; calculate projection circle center
		//onto ray and see if distance is > 0
		//Vector4d circleCtrProjection

		list<Vector4d> returnList;
		float a = inputRay->direction_.dot(inputRay->direction_);
		float b = 2 * inputRay->direction_.dot(inputRay->origin_ - centerPt_);
		float c = (inputRay->origin_ - centerPt_).dot(inputRay->origin_ - centerPt_) - radius_*radius_;
		float discriminant = pow(b,2) - 4*a*c;
		float quadraticFormResult = float(-1 * b - sqrt(b*b-4*a*c)) / float(2*a);

		Vector4d projectionPt = 
			((centerPt_.dot(inputRay->origin_)) / pow(centerPt_.norm(), 2)) * centerPt_;

		if (discriminant < 0) {
			return returnList;
		} 

		else if (discriminant == 0) {
			returnList.push_back(projectionPt);
		} else {
			Vector4d a = radius_;
			Vector4d b = projectionPt - centerPt_;
			b = b.cwiseAbs();
			Vector4d c = sqrt(pow(a,2) - b * b);
			Vector4d distRayOrigToProjPt = projectonPt - centerPt_;
			distRayOrigToProjPt = distRayOrigToProjPt.cwiseAbs();
			Vector4d distRayOrigToCenterPt = distRayOrigToProjPt - c;

			// finding first inersection point
			Vector4d intersectionPt1 = inputRay->origin_ + inputRay->direction_ * distRayOrigToCenterPt;

			// finding second intersection point
			Vector4d halfCircleSegment = projectionPt_ - intersectionPt1;
			Vector4d intersectionPt2 = intersectionPt1 + 2*halfCircleSegment;

			returnList.push_back(intersectionPt1);
			returnList.push_back(intersectionPt2);
		}
		return returnList;
	}

};

int main(int argc, char *argv[]) {

	Sphere testSphere (Vector3d w(0,0,0,0), 200.0);

	Ray testRay (Vector3d x(0,2,0,0), Vector3d y(0,-1,0,0), 0.5);

	return testSphere.sphereIntersection(&testRay);



	//return 0;
}
