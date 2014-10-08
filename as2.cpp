#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

int main(int argc, char *argv[]) {

typedef Matrix<double, 3, 1> Color;

class Ray () {
public:
	Vector4d origin_;
	Vector4d direction_;
	float t_;
}

class BRDF() {
	double specularCoeffecient;
	float ambientColor;
	Color ambientColor;

}



	return 0;
}
