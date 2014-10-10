#include <getopt.h>
#include <stdexcept>
#include <iostream>
#include "options.h"

Options programOptions;
const struct option Options::getoptOptions[] = {
	{"help", 0, NULL, OPTION_HELP},
	{"output", 1, NULL, OPTION_OUTPUT},
	{"width", 1, NULL, OPTION_WIDTH},
	{"height", 1, NULL, OPTION_HEIGHT},
	{"bdepth", 1, NULL, OPTION_BOUNCE_DEPTH},
	{"intersection-only", 0, NULL, OPTION_INTERSECTION_ONLY},
	{0}
};

bool Options::parseCommandLine(int argc, char *argv[]) {
	int opt;
	while ((opt = getopt_long(argc, argv, "w:h:o:", Options::getoptOptions, NULL)) != -1) {
		switch (opt) {
			case OPTION_OUTPUT:
				outputFilename_ = optarg;
				break;
			case OPTION_INTERSECTION_ONLY:
				intersectionOnly_ = true;
				break;
			case OPTION_WIDTH:
			case OPTION_HEIGHT:
			{
				int& dest = opt == OPTION_WIDTH ? renderWidth_ : renderHeight_;
				try {
					dest = std::stoi(optarg);
				} catch (std::logic_error& e) {
					std::cerr << "Error: Width and/or height is invalid." << std::endl;
					return false;
				}
				if (dest <= 0) {
					std::cerr << "Error: Width and/or height must be positive." << std::endl;
					return false;
				}
				break;
			}
			case OPTION_BOUNCE_DEPTH:
				try {
					bounceDepth_ = std::stoi(optarg);
				} catch (std::logic_error& e) {
					std::cerr << "Error: Bounce depth is invalid." << std::endl;
					return false;
				}
				if (bounceDepth_ < 0) {
					std::cerr << "Error: Bounce depth must be non-negative." << std::endl;
					return false;
				}
				break;
			case OPTION_HELP:
			case '?':
				printHelp(argv[0]);
				return false;
		}
	}
	while (optind < argc) {
		inputFilenames_.push_back(argv[optind++]);
	}
	if (inputFilenames_.empty()) {
		std::cerr << "Error: At least one input file must be specified." << std::endl;
		return false;
	}
	if (outputFilename_.empty()) {
		std::cerr << "Error: An output file must be specified." << std::endl;
		return false;
	}
	return true;
}

void Options::printHelp(const char* prog) {
	std::cerr << "Usage: " << prog << " [options] -o <output file> <input files>..." << std::endl;
}
