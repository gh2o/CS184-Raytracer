#pragma once
#include <string>
#include <vector>

class Options {
public:
	bool parseCommandLine(int argc, char *argv[]);
	void printHelp(const char *prog);
public:
	std::vector<std::string> inputFilenames_;
	std::string outputFilename_;
	int renderWidth_ = 500;
	int renderHeight_ = 500;
	int bounceDepth_ = 10;
	bool intersectionOnly_ = false;
public:
	static const struct option getoptOptions[];
	enum {
		OPTION_HELP,
		OPTION_OUTPUT = 'o',
		OPTION_WIDTH = 'w',
		OPTION_HEIGHT = 'h',
		OPTION_BOUNCE_DEPTH,
		OPTION_INTERSECTION_ONLY
	};
};

extern Options programOptions;
