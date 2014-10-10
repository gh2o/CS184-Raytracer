#pragma once

class Util {
public:
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
