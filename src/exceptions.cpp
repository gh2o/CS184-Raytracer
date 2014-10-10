#include "exceptions.h"

std::string ParseException::buildMessage(std::string msg, int lineno) {
	if (lineno > 0) {
		std::ostringstream s;
		s << "line " << lineno << ": " << msg;
		return s.str();
	} else {
		return msg;
	}
}
