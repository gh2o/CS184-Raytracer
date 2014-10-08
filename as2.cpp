#include <getopt.h>
#include <png.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <locale>
#include <memory>
#include <map>
#include <Eigen/Dense>

using namespace Eigen;

typedef Vector3d Color;

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
public:
	ParseException(std::string msg) :
		ParseException(msg, -1) {}
	ParseException(std::string msg, int lineno) :
		runtime_error(buildMessage(msg, lineno)) {}
private:
	std::string buildMessage(std::string msg, int lineno) {
		if (lineno > 0) {
			std::ostringstream s;
			s << "line " << lineno << ": " << msg;
			return s.str();
		} else {
			return msg;
		}
	}
};

class WriteException : public std::runtime_error {
public:
	using runtime_error::runtime_error;
};

class GlobalScene {
public:
	typedef Array<Vector3d, Dynamic, Dynamic, RowMajor> RasterImage;
public:
	void renderScene(RasterImage& output) {
		for (int i = 200; i < 600; i++)
			for (int j = 200; j < 400; j++)
				output(i,j) << 0,1,1;
	}
};

class RTParser {
private:
	enum LineTypeValue {
		LINE_TYPE_CAMERA,
		LINE_TYPE_SPHERE,
		LINE_TYPE_TRIANGLE,
		LINE_TYPE_POINT_LIGHT,
		LINE_TYPE_DIRECTIONAL_LIGHT,
		LINE_TYPE_AMBIENT_LIGHT,
		LINE_TYPE_MATERIAL,
		LINE_TYPE_TRANSFORM_TRANSLATE,
		LINE_TYPE_TRANSFORM_ROTATE,
		LINE_TYPE_TRANSFORM_SCALE,
		LINE_TYPE_TRANSFORM_IDENTITY
	};
	struct LineType {
		LineTypeValue type_;
		int pmin_, pmax_;
		LineType() = default;
		LineType(const LineType&) = default;
		LineType(LineTypeValue type, int pval) :
			LineType(type, pval, pval) {}
		LineType(LineTypeValue type, int pmin, int pmax):
			type_(type), pmin_(pmin), pmax_(pmax) {}
	};
	static const std::map<std::string, LineType> LINE_TYPES;
	static std::map<std::string, LineType> initializeLineTypes() {
		return {
			{"cam", {LINE_TYPE_CAMERA, 15}},
			{"sph", {LINE_TYPE_SPHERE, 4}},
			{"tri", {LINE_TYPE_TRIANGLE, 9}},
			{"ltp", {LINE_TYPE_POINT_LIGHT, 6, 7}},
			{"ltd", {LINE_TYPE_DIRECTIONAL_LIGHT, 6}},
			{"lta", {LINE_TYPE_AMBIENT_LIGHT, 3}},
			{"mat", {LINE_TYPE_MATERIAL, 13}},
			{"xft", {LINE_TYPE_TRANSFORM_TRANSLATE, 3}},
			{"xfr", {LINE_TYPE_TRANSFORM_ROTATE, 3}},
			{"xfs", {LINE_TYPE_TRANSFORM_SCALE, 3}},
			{"xfz", {LINE_TYPE_TRANSFORM_IDENTITY, 0}}
		};
	}
public:
	RTParser(GlobalScene& scene) : scene_(scene) {}
	void parseFile(std::string filename) {
		std::ifstream stream(filename);
		if (!stream)
			throw ParseException("file not found: " + filename);
		int lineno = 1;
		for (std::string line; std::getline(stream, line); lineno++) {
			std::istringstream ss(line);
			std::string stype = extractToken(ss, lineno);
			if (stype.empty()) {
				// blank line
				continue;
			}
			if (stype == "obj") {
				// obj needs special handling
				std::string filename = extractToken(ss, lineno);
				if (filename.empty())
					throw ParseException("obj requires a filename", lineno);
				// TODO
				abort();
				continue;
			}
			auto iter = LINE_TYPES.find(stype);
			if (iter == LINE_TYPES.end()) {
				std::cerr << "Warning: unknown line type " << stype << std::endl;
				continue;
			}
			LineType ltype = iter->second;
			std::unique_ptr<double[]> params(new double[ltype.pmax_]);
			std::fill(params.get(), params.get() + ltype.pmax_, 0.0);
			for (int i = 0; i < ltype.pmax_; i++) {
				std::string token = extractToken(ss, lineno);
				if (token.empty()) {
					if (i < ltype.pmin_) {
						std::ostringstream es;
						es << stype
							<< " requires "
							<< (ltype.pmin_ == ltype.pmax_ ? "" : "at least ")
							<< ltype.pmin_
							<< " parameters";
						throw ParseException(es.str(), lineno);
					} else {
						break;
					}
				}
				try {
					params[i] = std::stod(token);
				} catch (std::logic_error& e) {
					std::ostringstream es;
					es << "invalid number " << token;
					throw ParseException(es.str(), lineno);
				}
			}
			// TODO
			abort();
		}
	}
	std::string extractToken(std::istream& stream, int lineno) {
		// skip beginning spaces
		int c;
		while (true) {
			c = stream.peek();
			if (c == EOF)
				return std::string();
			if (!std::isspace(c))
				break;
			stream.get();
		}
		// read string value
		std::ostringstream s;
		if (stream.peek() == '"') {
			stream.get();
			while (true) {
				c = stream.get();
				if (c == EOF)
					throw ParseException("unclosed quotes", lineno);
				if (c == '"')
					break;
				s.put(c);
			}
		} else {
			while (true) {
				c = stream.get();
				if (c == EOF || std::isspace(c))
					break;
				s.put(c);
			}
		}
		// return it
		return s.str();
	}
private:
	GlobalScene& scene_;
};

class PNGWriter {
private:
	typedef Matrix<uint8_t, 3, 1> Vector3uc;
	typedef Array<Vector3uc, Dynamic, Dynamic, RowMajor> RGBImage;
public:
	PNGWriter(std::string filename) : filename_(filename) {}
	void writeImage(GlobalScene::RasterImage& image) {
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
	RGBImage convertImage(GlobalScene::RasterImage& image) {
		RGBImage rgb(image.rows(), image.cols());
		for (int i = 0; i < image.size(); i++)
			rgb(i) = (image(i).cwiseMin(1).cwiseMax(0) * 255.0).cast<uint8_t>();
		return rgb;
	}
private:
	std::string filename_;
};

const std::map<std::string, RTParser::LineType> RTParser::LINE_TYPES(RTParser::initializeLineTypes());

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
	// check that output is writable
	if (!std::ofstream(outputFilename)) {
		std::cerr << "Error: Output file is not writable." << std::endl;
		return 1;
	} else {
		std::remove(outputFilename.c_str());
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
	GlobalScene::RasterImage image(2000,3000);
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
