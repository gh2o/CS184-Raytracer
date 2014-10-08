#include <getopt.h>
#include <png.h>
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
	Color diffuseColor;
	Color specularColor;
	Color reflectiveColor;
} Material;

class ParseException : public std::runtime_error {
};

class GlobalScene {
public:
	typedef Matrix<Vector3d, Dynamic, Dynamic> RasterImage;
public:
	void renderScene(RasterImage& output) {
	}
};

class RTParser {
public:
	RTParser(GlobalScene& scene) : scene_(scene) {}
	void parseFile(std::string filename) {
	}
private:
	GlobalScene& scene_;
};

class PNGWriter {
public:
	PNGWriter(std::string filename) : filename_(filename) {}
	void writeImage(GlobalScene::RasterImage& image) {
	}
private:
	std::string filename_;
};

static struct option programOptions[] = {
	{"help", 0, NULL, 'h'},
	{"output", 1, NULL, 'o'}
};

static void printHelp(const char *prog) {
	std::cerr << "Usage: " << prog << " [options] -o <output file> <input files>..." << std::endl;
}

int main(int argc, char *argv[]) {
	// program options
	std::string outputFilename;
	// parse options
	int opt;
	while ((opt = getopt_long(argc, argv, "ho:", programOptions, NULL)) != -1) {
		switch (opt) {
			case 'o':
				outputFilename = optarg;
				std::cerr << "output: " << optarg << std::endl;
				break;
			case 'h':
			case '?':
				printHelp(argv[0]);
				return 1;
		}
	}
	if (optind >= argc) {
		std::cerr << "Error: At least one input file must be specified." << std::endl;
		return 1;
	}
	if (outputFilename.empty()) {
		std::cerr << "Error: An output file must be specified." << std::endl;
		return 1;
	}
	// read input
	GlobalScene scene;
	while (optind < argc) {
		RTParser parser(scene);
		try {
			parser.parseFile(argv[optind++]);
		} catch (const ParseException& e) {
			std::cerr << "Error: " << e.what() << std::endl;
			return 1;
		}
	}
	// render scene
	GlobalScene::RasterImage image;
	scene.renderScene(image);
	// write output
	PNGWriter writer(outputFilename);
	writer.writeImage(image);
	// done!
	return 0;
}
