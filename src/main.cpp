#include <getopt.h>
#include <png.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <locale>
#include <memory>
#include <map>
#include <vector>
#include <array>
#include <limits>
#include "common.h"
#include "exceptions.h"
#include "options.h"
#include "util.h"
#include "rtbase.h"
#include "lights.h"
#include "geometry.h"
#include "scene.h"
#include "parsers.h"

class PNGWriter {
private:
	typedef Matrix<uint8_t, 3, 1> Vector3uc;
	typedef Array<Vector3uc, Dynamic, Dynamic, RowMajor> RGBImage;
public:
	PNGWriter(std::string filename) : filename_(filename) {}
	void writeImage(Scene::RasterImage& image) {
		RGBImage rgb = convertImage(image);
		png_image png = {0};
		png.version = PNG_IMAGE_VERSION;
		png.width = image.cols();
		png.height = image.rows();
		png.format = PNG_FORMAT_RGB;
		if (!png_image_write_to_file(&png, filename_.c_str(), false,
		                             rgb.data(), image.cols() * 3, NULL))
			throw WriteException(std::string(png.message));
	}
	RGBImage convertImage(Scene::RasterImage& image) {
		RGBImage rgb(image.rows(), image.cols());
		for (int i = 0; i < image.size(); i++)
			rgb(i) = (image(i).cwiseMin(1).cwiseMax(0) * 255.0).cast<uint8_t>();
		return rgb;
	}
private:
	std::string filename_;
};

int main(int argc, char *argv[]) {
	// program options
	if (!programOptions.parseCommandLine(argc, argv))
		return 1;
	// check that output is writable
	const std::string& outputFilename = programOptions.outputFilename_;
	if (!std::ofstream(outputFilename)) {
		std::cerr << "Error: Output file is not writable." << std::endl;
		return 1;
	} else {
		std::remove(outputFilename.c_str());
	}
	// read input
	Scene scene;
	for (std::string inputFilename : programOptions.inputFilenames_) {
		RTParser parser(scene);
		try {
			parser.parseFile(inputFilename);
		} catch (const ParseException& e) {
			std::cerr << "Error: " << e.what() << std::endl;
			return 1;
		}
	}
	if (!scene.hasCamera()) {
		std::cerr << "Error: At least one camera must be specified." << std::endl;
		return 1;
	}
	// render scene
	Scene::RasterImage image(programOptions.renderHeight_,
	                               programOptions.renderWidth_);
	scene.renderScene(image);
	// write output
	try {
		PNGWriter writer(outputFilename);
		writer.writeImage(image);
	} catch (const WriteException& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return 1;
	}
	// done!
	return 0;
}
