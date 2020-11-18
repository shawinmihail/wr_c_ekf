#include "Utils.h"
#include <math.h>

const float UTILS_EPS = 1e-6f;

Vector4 quatFromEul(const Vector3& eul)
{
	 float roll = eul[0];
	 float pitch = eul[1];
	 float yaw = eul[2];

	 float cy = cos(yaw * 0.5f);
	 float sy = sin(yaw * 0.5f);
	 float cr = cos(roll * 0.5f);
	 float sr = sin(roll * 0.5f);
	 float cp = cos(pitch * 0.5f);
	 float sp = sin(pitch * 0.5f);

	 float q0 = cy * cr * cp + sy * sr * sp;
	 float q1 = cy * sr * cp - sy * cr * sp;
	 float q2 = cy * cr * sp + sy * sr * cp;
	 float q3 = sy * cr * cp - cy * sr * sp;

	 return Vector4(q0, q1, q2 ,q3);
}

Vector3 quatToEul(const Vector4& q)
{
    double sinr = 2.0f * (q[0] * q[1] + q[2] * q[3]);
    double cosr = 1.0f - 2.0f * (q[1] * q[1] + q[2] * q[2]);
    double roll = atan2(sinr, cosr);

    double sinp = 2.0 * (q[0] * q[2] - q[3] * q[1]);
    double pitch = 0.;
    if (abs(sinp) >= 1) {
        pitch = sinp / fabs (sinp) * 3.1415 / 2.0;
    }
    else {
        pitch = asin(sinp);
    }

    double siny = 2.0 * (q[0] * q[3] + q[1] * q[2]);
    double cosy = 1.0 - 2.0 * (q[2] * q[2] + q[3] * q[3]);
    double yaw = atan2(siny, cosy);
    return Vector3(roll, pitch, yaw);
}

Vector4 quatMultiply(const Vector4& q, const Vector4& r)
{
	Vector4 p;
	p[0] = r[0] * q[0] - r[1] * q[1] - r[2] * q[2] - r[3] * q[3];
	p[1] = r[0] * q[1] + r[1] * q[0] - r[2] * q[3] + r[3] * q[2];
	p[2] = r[0] * q[2] + r[1] * q[3] + r[2] * q[0] - r[3] * q[1];
	p[3] = r[0] * q[3] - r[1] * q[2] + r[2] * q[1] + r[3] * q[0];
	return p;
}

Vector4 quatInverse(const Vector4& q)
{
	Vector4 qDual (q[0], -q[1], -q[2], -q[3]);
	return qDual;
}

Vector3 quatRotate(const Vector4& q, const Vector3& v)
{
	Vector4 qv(0.0f, v[0], v[1], v[2]);

	Vector4 qDual = quatInverse(q);
	Vector4 qv1 = quatMultiply(qv, qDual);
	Vector4 qv2 = quatMultiply(q, qv1);

	return Vector3(qv2[1], qv2[2], qv2[3]);
}

Eigen::Matrix<float, 3, 3> quatToMatrix(const Vector4& q)
{
	float R11 = 1.f - 2.f * q(2) * q(2) - 2.f * q(3) * q(3);
	float R12 = 2.f * q(1) * q(2) - 2.f * q(3) * q(0);
	float R13 = 2.f * q(1) * q(3) + 2.f * q(2) * q(0);

	float R21 = 2.f * q(1) * q(2) + 2.f * q(3) * q(0);
	float R22 = 1.f - 2.f * q(1) * q(1) - 2.f * q(3) * q(3);
	float R23 = 2.f * q(2) * q(3) - 2.f * q(1) * q(0);

	float R31 = 2.f * q(1) * q(3) - 2.f * q(2) * q(0);
	float R32 = 2.f * q(2) * q(3) + 2.f * q(1) * q(0);
	float R33 = 1.f - 2.f * q(1) * q(1) - 2.f * q(2) * q(2);

	Eigen::Matrix<float, 3, 3> R;
	R << R11, R12, R13, R21, R22, R23, R31, R32, R33;
	return R;
}

Eigen::Matrix<float, 3, 3> crossOperator(const Vector3& v)
{
	Eigen::Matrix<float, 3, 3> R;
	R << 0.f, -v(2), v(1),
		 v(2), 0.f, -v(0),
		-v(1), v(0), 0.f;
	return R;
}

Vector4 quatBetweenVectors(const Vector3& v1, const Vector3& v2)
{
    double eps = 1e-6;
    if ((v1.norm() < eps) || (v2.norm() < eps)) 
    {
        return Vector4(1.f, 0.f, 0.f, 0.f);
    }
    
    Vector3 pin = v1.cross(v2);
    double w = sqrtf((v1.norm() * v1.norm()) * (v2.norm() * v2.norm())) + v1.dot(v2);
    
    Vector4 q(w, pin[0], pin[1], pin[2]);
    if (q.norm() < eps)
    {
        return Vector4 (0.f,0.f,0.f,1.f);
    }
    
    q.normalize();
    return q;
}


Eigen::Matrix<float, 3, 4> quatRotateLinearizationQ(const Vector4& q, const Vector3& v)
{
	// f = q * v * q_dual
	// res = df/dq

	/*
	rddot_q_j(q1, q2, q3, q4) =
	[2 * a1 * q1 - 2 * a2 * q4 + 2 * a3 * q3, 2 * a1 * q2 + 2 * a2 * q3 + 2 * a3 * q4, 2 * a2 * q2 - 2 * a1 * q3 + 2 * a3 * q1, 2 * a3 * q2 - 2 * a1 * q4 - 2 * a2 * q1]
	[2 * a2 * q1 + 2 * a1 * q4 - 2 * a3 * q2, 2 * a1 * q3 - 2 * a2 * q2 - 2 * a3 * q1, 2 * a1 * q2 + 2 * a2 * q3 + 2 * a3 * q4, 2 * a1 * q1 - 2 * a2 * q4 + 2 * a3 * q3]
	[2 * a2 * q2 - 2 * a1 * q3 + 2 * a3 * q1, 2 * a2 * q1 + 2 * a1 * q4 - 2 * a3 * q2, 2 * a2 * q4 - 2 * a1 * q1 - 2 * a3 * q3, 2 * a1 * q2 + 2 * a2 * q3 + 2 * a3 * q4]
	*/

	Eigen::Matrix<float, 3, 4> res;
	res << v(0) * q(0) - v(1) * q(3) + v(2) * q(2), v(0)* q(1) + (2) * v(1) * q(2) + v(2) * q(3), v(1)* q(1) - v(0) * q(2) + v(2) * q(0), v(2)* q(1) - v(0) * q(3) - v(1) * q(0),
		   v(1) * q(0) + v(0) * q(3) - v(2) * q(1), v(0)* q(2) - (1) * v(1) * q(1) - v(2) * q(0), v(0)* q(1) + v(1) * q(2) + v(2) * q(3), v(0)* q(0) - v(1) * q(3) + v(2) * q(2),
		   v(1) * q(1) - v(0) * q(2) + v(2) * q(0), v(1)* q(0) + (1) * v(0) * q(3) - v(2) * q(1), v(1)* q(3) - v(0) * q(0) - v(2) * q(2), v(0)* q(1) + v(1) * q(2) + v(2) * q(3);

	return 2.f*res;
}

Eigen::Matrix<float, 4, 4> poissonEqLinearizationQ(const Vector3& w)
{
	// f = 0.5 * q0 * qw0
	// res = df/dq

	/*
	qdot_q_j(q1, q2, q3, q4) =
	[    0, -w1/2, -w2/2, -w3/2]
	[ w1/2,     0,  w3/2, -w2/2]
	[ w2/2, -w3/2,     0,  w1/2]
	[ w3/2,  w2/2, -w1/2,     0]
	*/

	Eigen::Matrix<float, 4, 4> res;
	res << 0, -w(0) / 2.f, -w(1) / 2.f, -w(2) / 2.f,
		  w(0) / 2.f, 0, w(2) / 2.f, -w(1) / 2.f,
		  w(1) / 2.f, -w(2) / 2.f, 0, w(0) / 2.f,
		  w(2) / 2.f, w(1) / 2.f, -w(0) / 2.f, 0;

	return res;
}

Eigen::Matrix<float, 4, 3> poissonEqLinearizationW(const Vector4& q, const Vector3& w)
{
	// f = 0.5 * q0 * qw0
	// res = df/dq

	/*
	qdot_w_j(w1, w2, w3) =

	[ -q2/2, -q3/2, -q4/2]
	[  q1/2, -q4/2,  q3/2]
	[  q4/2,  q1/2, -q2/2]
	[ -q3/2,  q2/2,  q1/2]
	*/

	Eigen::Matrix<float, 4, 3> res;
	res << -q(1) / 2.f, -q(2) / 2.f, -q(3) / 2.f,
		    q(0) / 2.f, -q(3) / 2.f, q(2) / 2.f,
			q(3) / 2.f, q(0) / 2.f, -q(1) / 2.f,
			-q(2) / 2.f, q(1) / 2.f, q(0) / 2.f;

	return res;
}

Eigen::Matrix<float, 1, 3> normVect3Linearization(const Vector3& v)
{
	Eigen::Matrix<float, 1, 3> res;
	float vn = v.norm();
	res << v(0) / vn, v(1) / vn, v(2) / vn;
	return res;
}

