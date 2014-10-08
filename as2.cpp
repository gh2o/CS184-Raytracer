#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

typedef Matrix<double, 3, 1> Color;

struct {
	Vector4d origin_;
	Vector4d direction_;
	float t_;
} Ray;

struct {
	Color ambientColor;
	Color ambientColor;
	Color specularColor;
	Color reflectiveColor;
} Material;

int main(int argc, char *argv[]) {


	return 0;
}
