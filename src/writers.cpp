#include <png.h>
#include "writers.h"

PNGWriter::RGBImage PNGWriter::convertToRGBImage(Scene::RasterImage& image) {
	RGBImage rgb(image.rows(), image.cols());
	for (int i = 0; i < image.size(); i++)
		rgb(i) = (image(i).cwiseMin(1).cwiseMax(0) * 255.0).cast<uint8_t>();
	return rgb;
}

void PNGWriter::writeImage(Scene::RasterImage& image) {
	RGBImage rgb = convertToRGBImage(image);
	png_image png = {0};
	png.version = PNG_IMAGE_VERSION;
	png.width = image.cols();
	png.height = image.rows();
	png.format = PNG_FORMAT_RGB;
	if (!png_image_write_to_file(&png, filename_.c_str(), false,
								 rgb.data(), image.cols() * 3, NULL))
		throw WriteException(std::string(png.message));
}
