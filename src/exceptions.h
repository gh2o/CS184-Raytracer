#pragma once
#include <stdexcept>
#include <iostream>
#include <sstream>

class ParseException : public std::runtime_error {
public:
	ParseException(std::string msg) :
		ParseException(msg, -1) {}
	ParseException(std::string msg, int lineno) :
		runtime_error(buildMessage(msg, lineno)) {}
	static void showWarning(std::string msg) {
		showWarning(msg, -1);
	}
	static void showWarning(std::string msg, int lineno) {
		std::cerr << "Warning: " << buildMessage(msg, lineno) << std::endl;
	}
private:
	static std::string buildMessage(std::string msg, int lineno);
};

class MathException : public std::runtime_error {
public:
	using runtime_error::runtime_error;
};

class WriteException : public std::runtime_error {
public:
	using runtime_error::runtime_error;
};
