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

float shortestRotation(float from, float to)
{
    float pi = 3.1415;
    float da = fmod((to - from + pi), 2*pi) - pi;
    if (da < -pi)
    {
        da = da + 2*pi;
    }
    return da;
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

Vector3 eulKinematics(const Vector3& e, const Vector3& w)
{
    // Need to check on theta +- pi/2
    float psi = e[2];
    float theta = e[1];
    float phi = e[0];
    Eigen::Matrix<float, 3, 3> A;
    A <<   1,   sin(phi)*tan(theta),   cos(phi)*tan(theta),
           0,   cos(phi),              -sin(phi),
           0,   sin(phi)/cos(theta),   cos(phi)/cos(theta);

    Vector3 eul_dot = A*w;
    return eul_dot;
}

Eigen::Matrix<float, 3, 3> eulToMatrix(const Vector3& e)
{
    float sa = sin(e[2]);
    float ca = cos(e[2]);
    float sb = sin(e[1]);
    float cb = cos(e[1]);
    float sg = sin(e[0]);
    float cg = cos(e[0]);
    
    float R11 = ca * cb;
    float R21 = sa * cb;
    float R31 = -sb;

    float R12 = ca * sb * sg  - sa * cg;
    float R22 = sa * sb * sg + ca * cg;
    float R32 = cb * sg;

    float R13 = ca * sb * cg + sa *sg;
    float R23 = sa * sb * cg - ca * sg;
    float R33 = cb * cg;

    Eigen::Matrix<float, 3, 3> R;
    R <<   R11,   R12,   R13,
           R21,   R22,   R23,
           R31,   R32,   R33;

    return R;
}

Eigen::Matrix<float, 3, 3> eulRotLinearization(const Vector3& e, const Vector3& w)
{
    
    float se0 = sin(e[0]);
    float se1 = sin(e[1]);
    float se2 = sin(e[2]);
    float ce0 = cos(e[0]);
    float ce1 = cos(e[1]);
    float ce2 = cos(e[2]);
    
    float A11 = w[1]*(se0*se2 + ce0*ce2*se1) + w[2]*(ce0*se2 - ce2*se0*se1);
    float A12 = w[2]*ce0*ce1*ce2 - w[0]*ce2*se1 + w[1]*ce1*ce2*se0;
    float A13 = w[2]*(ce2*se0 - ce0*se1*se2) - w[1]*(ce0*ce2 + se0*se1*se2) - w[0]*ce1*se2;

    float A21 = -w[1]*(ce2*se0 - ce0*se1*se2) - w[2]*(ce0*ce2 + se0*se1*se2);
    float A22 = w[2]*ce0*ce1*se2 - w[0]*se1*se2 + w[1]*ce1*se0*se2;
    float A23 = w[2]*(se0*se2 + ce0*ce2*se1) - w[1]*(ce0*se2 - ce2*se0*se1) + w[0]*ce1*ce2;

    float A31 = w[1]*ce0*ce1 - w[2]*ce1*se0;
    float A32 = -w[0]*ce1 - w[2]*ce0*se1 -w[1]*se0*se1;
    float A33 = 0;
    
    Eigen::Matrix<float, 3, 3> A;
    A <<   A11,   A12,   A13,
           A21,   A22,   A23,
           A31,   A32,   A33;

    return A;
}

Eigen::Matrix<float, 3, 3> eulKinLinearization(const Vector3& e, const Vector3& w)
{
    float c_phi = cos(e[2]);
    float s_phi = sin(e[2]);
    float t_theta = tan(e[1]);
    float c_theta = cos(e[1]);
    float s_theta = sin(e[1]);

    float A11 = w[1]*c_phi*t_theta - w[2]*s_phi*t_theta;
    float A12 = w[2]*c_phi*(t_theta*t_theta + 1) + w[1]*s_phi*(t_theta*t_theta + 1);
    float A13 = 0;

    float A21 = - w[2]*c_phi - w[1]*s_phi;
    float A22 = 0;
    float A23 = 0;

    float A31 = (w[1]*c_phi)/c_theta - (w[2]*s_phi)/c_theta;
    float A32 = (w[2]*c_phi*s_theta)/c_theta*c_theta + (w[1]*s_phi*s_theta)/c_theta*c_theta;
    float A33 = 0;

    Eigen::Matrix<float, 3, 3> A;
    A <<   A11,   A12,   A13,
           A21,   A22,   A23,
           A31,   A32,   A33;

    return A;
}

