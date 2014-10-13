#pragma once
#include <libgen.h>

class Util {
public:
	template<typename Func>
	static std::string genname(const std::string& path, Func func) {
		auto length = path.length();
		char* buffer = new char[length + 1];
		path.copy(buffer, length);
		buffer[length] = '\0';
		std::string result = func(buffer);
		delete[] buffer;
		return result;
	}
	static std::string basename(const std::string& path) {
		return genname(path, ::basename);
	}
	static std::string dirname(const std::string& path) {
		return genname(path, ::dirname);
	}
	static Vector4d cross(const Vector4d& a, const Vector4d& b) {
		Vector4d c = Vector4d::Zero();
		c.head<3>() = a.head<3>().cross(b.head<3>());
		return c;
	}
	/* NOTE: vec?dFrom?d functions should only be used with
	 *       directional (non-homogeneous) vectors. */
	static Vector4d vec4dFrom3d(const Vector3d& f) {
		Vector4d t = Vector4d::Zero();
		t.head<3>() = f;
		return t;
	}
	static Vector3d vec3dFrom4d(const Vector4d& f) {
		return f.head<3>();
	}
};
