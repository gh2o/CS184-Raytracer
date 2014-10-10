#pragma once
#include "common.h"
#include "scene.h"

class PNGWriter {
private:
	typedef Matrix<uint8_t, 3, 1> Vector3uc;
	typedef Array<Vector3uc, Dynamic, Dynamic, RowMajor> RGBImage;
public:
	PNGWriter(std::string filename) : filename_(filename) {}
	void writeImage(Scene::RasterImage& image);
private:
	static RGBImage convertToRGBImage(Scene::RasterImage& image);
private:
	std::string filename_;
};
