#include <signal.h>
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include "options.h"
#include "scene.h"
#include "parsers.h"
#include "writers.h"

static volatile bool progressSignaled = false;

static void updateProgress(int complete, int total) {
	if (complete != total && !progressSignaled)
		return;
	std::putchar('\r');
	std::printf("Rendering scene (%d/%d) (%.1f%%) ...",
			complete, total, 100.0 * complete / total);
	std::fflush(stdout);
	if (complete == total)
		std::putchar('\n');
	progressSignaled = false;
}

static void alarmProgress(int sig) {
	progressSignaled = true;
}

static void setupAlarm(bool enable) {
	if (enable)
		signal(SIGALRM, alarmProgress);
	suseconds_t period = enable ? 1000000 / 4 : 0;
	struct itimerval itv = {{0}};
	itv.it_value.tv_usec = period;
	itv.it_interval.tv_usec = period;
	setitimer(ITIMER_REAL, &itv, nullptr);
	if (!enable)
		signal(SIGALRM, SIG_DFL);
}

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
	signal(SIGALRM, alarmProgress);
	setupAlarm(true);
	Scene::RasterImage image(programOptions.renderHeight_,
	                               programOptions.renderWidth_);
	scene.renderScene(image, updateProgress);
	setupAlarm(false);
	signal(SIGALRM, SIG_DFL);
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
