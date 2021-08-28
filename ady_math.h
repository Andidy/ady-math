#ifndef ADY_MATH_H
#define ADY_MATH_H

#include <stdint.h>
#include <math.h>

inline int32_t intmax(int32_t a, int32_t b) {
	return (a > b) ? a : b;
}

inline int32_t intmin(int32_t a, int32_t b) {
	return (a < b) ? a : b;
}

inline float floatmax(float a, float b) {
	return (a > b) ? a : b;
}

inline float floatmin(float a, float b) {
	return (a < b) ? a : b;
}

const float PI = 3.14159265359f;

inline float DegToRad(float angle) {
	return (angle * (PI / 180.0f));
}

// Vector 2 ----------------------------------------------------------------

union vec2 {
	struct { float x; float y; };
	struct { float u; float v; };
	float data[2];
};

inline vec2 Vec2(float x, float y) {
	return { x, y };
}

union ivec2 {
	struct { int32_t x; int32_t y; };
	struct { int32_t u; int32_t v; };
	int32_t data[2];
};

inline ivec2 IVec2(int32_t x, int32_t y) {
	return { x, y };
}

inline int32_t Vec2To1DIndex(vec2 v, int32_t width) {
	return (int32_t)v.x + (int32_t)v.y * width;
}

inline bool Vec2Equals(vec2 v1, vec2 v2) {
	return (bool)((v1.x == v2.x) && (v1.y == v2.y));
}

// Vector 3 ----------------------------------------------------------------

union vec3 {
	struct { float x; float y; float z; };
	struct { float u; float v; float w; };
	struct { float r; float g; float b; };
	float data[3];
};

inline vec3 Vec3(float x, float y, float z) {
	return { x, y, z };
}

inline vec3 ZeroVec() {
	return { 0.0f, 0.0f, 0.0f };
}

inline vec3 OneVec() {
	return { 1.0f, 1.0f, 1.0f };
}

inline vec3 UpVec() {
	return { 0.0f, 1.0f, 0.0f };
}

inline bool VecEquals(vec3 v1, vec3 v2) {
	return (bool)((v1.x == v2.x) && (v1.y == v2.y) && (v1.z == v2.z));
}

inline vec3 NegVec(vec3 v1) {
	return { -1.0f * v1.x, -1.0f * v1.y, -1.0f * v1.z };
}

inline vec3 AddVec(vec3 v1, vec3 v2) {
	return { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
}

inline vec3 SubVec(vec3 v1, vec3 v2) {
	return { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
}

inline vec3 ScaleVec(vec3 v1, float f) {
	return { v1.x * f, v1.y * f, v1.z * f };
}

inline vec3 MulVec(vec3 v1, vec3 v2) {
	return { v1.x * v2.x, v1.y * v2.y, v1.z * v2.z };
}

inline float Dot(vec3 v1, vec3 v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline vec3 Cross(vec3 v1, vec3 v2) {
	return { v1.y * v2.z - v1.z * v2.y,
			 v1.z * v2.x - v1.x * v2.z,
			 v1.x * v2.y - v1.y * v2.x };
}

inline float Distance(vec3 v1, vec3 v2) {
	vec3 dist = SubVec(v2, v1);
	return sqrtf(Dot(dist, dist));
}

// TODO: Replace sqrtf with custom implementation?
inline float VecLen(vec3 v) {
	return sqrtf(Dot(v, v));
}

inline vec3 NormVec(vec3 v) {
	float len = VecLen(v);
	len = (len == 0.0f) ? 1.0f : 1.0f / len;
	return ScaleVec(v, len);
}

// Vector 4 ----------------------------------------------------------------

union vec4 {
	struct { float x; float y; float z; float w; };
	struct { float r; float g; float b; float a; };
	float data[4];
};

inline vec4 Vec4(float x, float y, float z, float w) {
	return { x, y, z, w };
}

// Matrix ------------------------------------------------------------------
// Originally styled after raysan5's raymath, and HandmadeMath

// Matrix functions and storage are row-major and right-handed coordinate system.
// When using matrices with vectors, consider the vectors as column vectors.
// This means that when multiplying matrices by vectors that it takes the form:
// M3 * M2 * M1 * v

// reminder of the memory layout
//
// [0][0] [0][1] [0][2] [0][3]
// [1][0] [1][1] [1][2] [1][3]
// [2][0] [2][1] [2][2] [2][3]
// [3][0] [3][1] [3][2] [3][3]
//
struct mat4 {
	float data[4][4];
};

inline vec4 MulMatVec(mat4 m, vec4 v) {
	vec4 result = Vec4(0.0f, 0.0f, 0.0f, 0.0f);

	for (int row = 0; row < 4; row++) {
		float sum = 0.0f;
		for (int col = 0; col < 4; col++) {
			sum += m.data[row][col] * v.data[col];
		}
		result.data[row] = sum;
	}

	return result;
}

inline mat4 MulMat(mat4 m2, mat4 m1) {
	mat4 result = {};
	for (int row = 0; row < 4; row++) {
		for (int col = 0; col < 4; col++) {
			float sum = 0;
			for (int i = 0; i < 4; i++) {
				sum += m2.data[row][i] * m1.data[i][col];
			}
			result.data[row][col] = sum;
		}
	}
	return result;
}

inline mat4 DiagonalMat(float diag) {
	mat4 result = {};

	result.data[0][0] = diag;
	result.data[1][1] = diag;
	result.data[2][2] = diag;
	result.data[3][3] = diag;

	return result;
}

inline mat4 TransposeMat(mat4 m) {
	mat4 result = {};

	result.data[0][0] = m.data[0][0];
	result.data[0][1] = m.data[1][0];
	result.data[0][2] = m.data[2][0];
	result.data[0][3] = m.data[3][0];

	result.data[1][0] = m.data[0][1];
	result.data[1][1] = m.data[1][1];
	result.data[1][2] = m.data[2][1];
	result.data[1][3] = m.data[3][1];

	result.data[2][0] = m.data[0][2];
	result.data[2][1] = m.data[1][2];
	result.data[2][2] = m.data[2][2];
	result.data[2][3] = m.data[3][2];

	result.data[3][0] = m.data[0][3];
	result.data[3][1] = m.data[1][3];
	result.data[3][2] = m.data[2][3];
	result.data[3][3] = m.data[3][3];

	return result;
}

inline mat4 ScaleMat(vec3 v) {
	mat4 result = DiagonalMat(1.0f);

	result.data[0][0] = v.x;
	result.data[1][1] = v.y;
	result.data[2][2] = v.z;

	return result;
}

inline mat4 TranslateMat(vec3 v) {
	mat4 result = DiagonalMat(1.0f);

	result.data[0][3] = v.x;
	result.data[1][3] = v.y;
	result.data[2][3] = v.z;

	return result;
}

// Note: angle is in degrees
inline mat4 RotateMatX(float angle) {
	mat4 result = DiagonalMat(1.0f);

	result.data[0][0] = 1.0f;

	float cos = cosf(DegToRad(angle));
	float sin = sinf(DegToRad(angle));

	result.data[1][1] = cos;
	result.data[1][2] = -sin;
	result.data[2][1] = sin;
	result.data[2][2] = cos;

	return result;
}

// Note: angle is in degrees
inline mat4 RotateMatY(float angle) {
	mat4 result = DiagonalMat(1.0f);

	float cos = cosf(DegToRad(angle));
	float sin = sinf(DegToRad(angle));

	result.data[0][0] = cos;
	result.data[0][2] = sin;
	result.data[1][1] = 1.0f;
	result.data[2][0] = -sin;
	result.data[2][2] = cos;

	return result;
}

// Note: angle is in degrees
inline mat4 RotateMatZ(float angle) {
	mat4 result = DiagonalMat(1.0f);

	float cos = cosf(DegToRad(angle));
	float sin = sinf(DegToRad(angle));

	result.data[0][0] = cos;
	result.data[0][1] = -sin;
	result.data[1][0] = sin;
	result.data[1][1] = cos;
	result.data[2][2] = 1.0f;

	return result;
}

// Note: angle is in degrees, axis will be normalized internally
inline mat4 RotateMat(float angle, vec3 axis) {
	mat4 result = DiagonalMat(1.0f);
	axis = NormVec(axis);

	float sin = sinf(DegToRad(angle));
	float cos = cosf(DegToRad(angle));
	float inv_cos = 1.0f - cos;

	result.data[0][0] = (axis.x * axis.x * inv_cos) + cos;
	result.data[0][1] = (axis.x * axis.y * inv_cos) - (axis.z * sin);
	result.data[0][2] = (axis.x * axis.z * inv_cos) + (axis.y * sin);

	result.data[1][0] = (axis.y * axis.x * inv_cos) + (axis.z * sin);
	result.data[1][1] = (axis.y * axis.y * inv_cos) + cos;
	result.data[1][2] = (axis.y * axis.z * inv_cos) - (axis.x * sin);

	result.data[2][0] = (axis.z * axis.x * inv_cos) - (axis.y * sin);
	result.data[2][1] = (axis.z * axis.y * inv_cos) + (axis.x * sin);
	result.data[2][2] = (axis.z * axis.z * inv_cos) + cos;

	return result;
}

inline mat4 LookAtMat(vec3 eye, vec3 target, vec3 up) {
	mat4 result = {};

	vec3 n = NormVec(SubVec(eye, target));
	vec3 u = NormVec(Cross(up, n));
	vec3 v = Cross(n, u);

	vec3 neg_eye = NegVec(eye);

	result.data[0][0] = u.x;
	result.data[0][1] = u.y;
	result.data[0][2] = u.z;
	result.data[0][3] = Dot(neg_eye, u);

	result.data[1][0] = v.x;
	result.data[1][1] = v.y;
	result.data[1][2] = v.z;
	result.data[1][3] = Dot(neg_eye, v);

	result.data[2][0] = n.x;
	result.data[2][1] = n.y;
	result.data[2][2] = n.z;
	result.data[2][3] = Dot(neg_eye, n);

	result.data[3][0] = 0.0f;
	result.data[3][1] = 0.0f;
	result.data[3][2] = 0.0f;
	result.data[3][3] = 1.0f;

	return result;
}

inline mat4 InverseLookAtMat(vec3 eye, vec3 target, vec3 up) {
	mat4 result = {};

	vec3 n = NormVec(SubVec(eye, target));
	vec3 u = NormVec(Cross(up, n));
	vec3 v = Cross(n, u);

	//vec3 neg_eye = NegVec(eye);

	result.data[0][0] = u.x;
	result.data[0][1] = v.x;
	result.data[0][2] = n.x;
	result.data[0][3] = eye.x;

	result.data[1][0] = u.y;
	result.data[1][1] = v.y;
	result.data[1][2] = n.y;
	result.data[1][3] = eye.y;

	result.data[2][0] = u.z;
	result.data[2][1] = v.z;
	result.data[2][2] = n.z;
	result.data[2][3] = eye.z;

	result.data[3][0] = 0.0f;
	result.data[3][1] = 0.0f;
	result.data[3][2] = 0.0f;
	result.data[3][3] = 1.0f;

	return result;
}

// Note: fov asks for an angle in degrees
inline mat4 PerspectiveMat(float fov, float aspect_ratio, float znear, float zfar) {
	mat4 result = {};

	float cotan = 1.0f / tanf(DegToRad(0.5f * fov));

	result.data[0][0] = cotan / aspect_ratio;
	result.data[1][1] = cotan;
	result.data[2][2] = -(zfar) / (zfar - znear);
	result.data[2][3] = -(znear * zfar) / (zfar - znear);
	result.data[3][2] = -1.0f;
	result.data[3][3] = 0.0f;

	return result;
}

inline mat4 OrthographicMat(float left, float right, float bot, float top, float znear, float zfar) {
	mat4 result = {};

	result.data[0][0] = 2.0f / (right - left);
	result.data[1][1] = 2.0f / (top - bot);
	result.data[2][2] = 2.0f / (znear - zfar);
	result.data[3][3] = 1.0f;

	result.data[3][0] = (left + right) / (left - right);
	result.data[3][1] = (bot + top) / (bot - top);
	result.data[3][2] = (zfar + znear) / (znear - zfar);

	return result;
}

inline mat4 ViewportMat(float width, float height) {
	mat4 result = {};

	result.data[0][0] = width / 2.0f;
	result.data[0][3] = width / 2.0f;
	result.data[1][1] = -height / 2.0f;
	result.data[1][3] = height / 2.0f;
	result.data[2][2] = 1.0f;
	result.data[2][3] = 0.0f;
	result.data[3][3] = 1.0f;

	return result;
}

inline mat4 ViewportPartialMat(float width, float height, float minX, float minY, float minZ, float maxZ) {
	mat4 result = {};

	result.data[0][0] = width / 2.0f;
	result.data[0][3] = minX + width / 2.0f;
	result.data[1][1] = -height / 2.0f;
	result.data[1][3] = minY + height / 2.0f;
	result.data[2][2] = maxZ - minZ;
	result.data[2][3] = minZ;
	result.data[3][3] = 1.0f;

	return result;
}

#endif