#include "EkfWgsMath.h"
#include <math.h>

const double UTILS_EPS = 1e-6f;

Eigen::Vector4d quatFromEul(const Eigen::Vector3d& eul)
{
	 double roll = eul[0];
	 double pitch = eul[1];
	 double yaw = eul[2];

	 double cy = cos(yaw * 0.5f);
	 double sy = sin(yaw * 0.5f);
	 double cr = cos(roll * 0.5f);
	 double sr = sin(roll * 0.5f);
	 double cp = cos(pitch * 0.5f);
	 double sp = sin(pitch * 0.5f);

	 double q0 = cy * cr * cp + sy * sr * sp;
	 double q1 = cy * sr * cp - sy * cr * sp;
	 double q2 = cy * cr * sp + sy * sr * cp;
	 double q3 = sy * cr * cp - cy * sr * sp;

	 return Eigen::Vector4d(q0, q1, q2 ,q3);
}

Eigen::Vector3d quatToEul(const Eigen::Vector4d& q)
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
    return Eigen::Vector3d(roll, pitch, yaw);
}

Eigen::Vector4d quatMultiply(const Eigen::Vector4d& q, const Eigen::Vector4d& r)
{
	Eigen::Vector4d p;
	p[0] = r[0] * q[0] - r[1] * q[1] - r[2] * q[2] - r[3] * q[3];
	p[1] = r[0] * q[1] + r[1] * q[0] - r[2] * q[3] + r[3] * q[2];
	p[2] = r[0] * q[2] + r[1] * q[3] + r[2] * q[0] - r[3] * q[1];
	p[3] = r[0] * q[3] - r[1] * q[2] + r[2] * q[1] + r[3] * q[0];
	return p;
}

Eigen::Vector4d quatInverse(const Eigen::Vector4d& q)
{
	Eigen::Vector4d qDual (q[0], -q[1], -q[2], -q[3]);
	return qDual;
}

Eigen::Vector3d quatRotate(const Eigen::Vector4d& q, const Eigen::Vector3d& v)
{
	Eigen::Vector4d qv(0.0f, v[0], v[1], v[2]);

	Eigen::Vector4d qDual = quatInverse(q);
	Eigen::Vector4d qv1 = quatMultiply(qv, qDual);
	Eigen::Vector4d qv2 = quatMultiply(q, qv1);

	return Eigen::Vector3d(qv2[1], qv2[2], qv2[3]);
}

Eigen::Matrix<double, 3, 3> quatToMatrix(const Eigen::Vector4d& q)
{
	double R11 = 1 - 2 * q(2) * q(2) - 2 * q(3) * q(3);
	double R12 = 2 * q(1) * q(2) - 2 * q(3) * q(0);
	double R13 = 2 * q(1) * q(3) + 2 * q(2) * q(0);

	double R21 = 2 * q(1) * q(2) + 2 * q(3) * q(0);
	double R22 = 1 - 2 * q(1) * q(1) - 2 * q(3) * q(3);
	double R23 = 2 * q(2) * q(3) - 2 * q(1) * q(0);

	double R31 = 2 * q(1) * q(3) - 2 * q(2) * q(0);
	double R32 = 2 * q(2) * q(3) + 2 * q(1) * q(0);
	double R33 = 1 - 2 * q(1) * q(1) - 2 * q(2) * q(2);

	Eigen::Matrix<double, 3, 3> R;
	R << R11, R12, R13, R21, R22, R23, R31, R32, R33;
	return R;
}

Eigen::Matrix<double, 3, 3> crossOperator(const Eigen::Vector3d& v)
{
	Eigen::Matrix<double, 3, 3> R;
	R << 0, -v(2), v(1),
		 v(2), 0, -v(0),
		-v(1), v(0), 0;
	return R;
}

Eigen::Vector4d quatBetweenVectors(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2)
{
    double eps = 1e-6;
    if ((v1.norm() < eps) || (v2.norm() < eps)) 
    {
        return Eigen::Vector4d(1, 0, 0, 0);
    }
    
    Eigen::Vector3d pin = v1.cross(v2);
    double w = sqrtf((v1.norm() * v1.norm()) * (v2.norm() * v2.norm())) + v1.dot(v2);
    
    Eigen::Vector4d q(w, pin[0], pin[1], pin[2]);
    if (q.norm() < eps)
    {
        return Eigen::Vector4d (0,0,0,1);
    }
    
    q.normalize();
    return q;
}

double shortestRotation(double from, double to)
{
    double pi = 3.1415;
    double da = fmod((to - from + pi), 2*pi) - pi;
    if (da < -pi)
    {
        da = da + 2*pi;
    }
    return da;
}


Eigen::Matrix<double, 3, 4> quatRotateLinearizationQ(const Eigen::Vector4d& q, const Eigen::Vector3d& v)
{
	// f = q * v * q_dual
	// res = df/dq

	/*
	rddot_q_j(q1, q2, q3, q4) =
	[2 * a1 * q1 - 2 * a2 * q4 + 2 * a3 * q3, 2 * a1 * q2 + 2 * a2 * q3 + 2 * a3 * q4, 2 * a2 * q2 - 2 * a1 * q3 + 2 * a3 * q1, 2 * a3 * q2 - 2 * a1 * q4 - 2 * a2 * q1]
	[2 * a2 * q1 + 2 * a1 * q4 - 2 * a3 * q2, 2 * a1 * q3 - 2 * a2 * q2 - 2 * a3 * q1, 2 * a1 * q2 + 2 * a2 * q3 + 2 * a3 * q4, 2 * a1 * q1 - 2 * a2 * q4 + 2 * a3 * q3]
	[2 * a2 * q2 - 2 * a1 * q3 + 2 * a3 * q1, 2 * a2 * q1 + 2 * a1 * q4 - 2 * a3 * q2, 2 * a2 * q4 - 2 * a1 * q1 - 2 * a3 * q3, 2 * a1 * q2 + 2 * a2 * q3 + 2 * a3 * q4]
	*/

	Eigen::Matrix<double, 3, 4> res;
	res << v(0) * q(0) - v(1) * q(3) + v(2) * q(2), v(0)* q(1) + (2) * v(1) * q(2) + v(2) * q(3), v(1)* q(1) - v(0) * q(2) + v(2) * q(0), v(2)* q(1) - v(0) * q(3) - v(1) * q(0),
		   v(1) * q(0) + v(0) * q(3) - v(2) * q(1), v(0)* q(2) - (1) * v(1) * q(1) - v(2) * q(0), v(0)* q(1) + v(1) * q(2) + v(2) * q(3), v(0)* q(0) - v(1) * q(3) + v(2) * q(2),
		   v(1) * q(1) - v(0) * q(2) + v(2) * q(0), v(1)* q(0) + (1) * v(0) * q(3) - v(2) * q(1), v(1)* q(3) - v(0) * q(0) - v(2) * q(2), v(0)* q(1) + v(1) * q(2) + v(2) * q(3);

	return 2*res;
}

Eigen::Matrix<double, 4, 4> poissonEqLinearizationQ(const Eigen::Vector3d& w)
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

	Eigen::Matrix<double, 4, 4> res;
	res << 0, -w(0) / 2, -w(1) / 2, -w(2) / 2,
		  w(0) / 2, 0, w(2) / 2, -w(1) / 2,
		  w(1) / 2, -w(2) / 2, 0, w(0) / 2,
		  w(2) / 2, w(1) / 2, -w(0) / 2, 0;

	return res;
}

Eigen::Matrix<double, 4, 3> poissonEqLinearizationW(const Eigen::Vector4d& q, const Eigen::Vector3d& w)
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

	Eigen::Matrix<double, 4, 3> res;
	res << -q(1) / 2, -q(2) / 2, -q(3) / 2,
		    q(0) / 2, -q(3) / 2, q(2) / 2,
			q(3) / 2, q(0) / 2, -q(1) / 2,
			-q(2) / 2, q(1) / 2, q(0) / 2;

	return res;
}

Eigen::Matrix<double, 1, 3> normVect3Linearization(const Eigen::Vector3d& v)
{
	Eigen::Matrix<double, 1, 3> res;
	double vn = v.norm();
	res << v(0) / vn, v(1) / vn, v(2) / vn;
	return res;
}

Eigen::Vector3d eulKinematics(const Eigen::Vector3d& e, const Eigen::Vector3d& w)
{
    // Need to check on theta +- pi/2
    double psi = e[2];
    double theta = e[1];
    double phi = e[0];
    Eigen::Matrix<double, 3, 3> A;
    A <<   1,   sin(phi)*tan(theta),   cos(phi)*tan(theta),
           0,   cos(phi),              -sin(phi),
           0,   sin(phi)/cos(theta),   cos(phi)/cos(theta);

    Eigen::Vector3d eul_dot = A*w;
    return eul_dot;
}

Eigen::Matrix<double, 3, 3> eulToMatrix(const Eigen::Vector3d& e)
{
    double sa = sin(e[2]);
    double ca = cos(e[2]);
    double sb = sin(e[1]);
    double cb = cos(e[1]);
    double sg = sin(e[0]);
    double cg = cos(e[0]);
    
    double R11 = ca * cb;
    double R21 = sa * cb;
    double R31 = -sb;

    double R12 = ca * sb * sg  - sa * cg;
    double R22 = sa * sb * sg + ca * cg;
    double R32 = cb * sg;

    double R13 = ca * sb * cg + sa *sg;
    double R23 = sa * sb * cg - ca * sg;
    double R33 = cb * cg;

    Eigen::Matrix<double, 3, 3> R;
    R <<   R11,   R12,   R13,
           R21,   R22,   R23,
           R31,   R32,   R33;

    return R;
}

Eigen::Matrix<double, 3, 3> eulRotLinearization(const Eigen::Vector3d& e, const Eigen::Vector3d& w)
{
    
    double se0 = sin(e[0]);
    double se1 = sin(e[1]);
    double se2 = sin(e[2]);
    double ce0 = cos(e[0]);
    double ce1 = cos(e[1]);
    double ce2 = cos(e[2]);
    
    double A11 = w[1]*(se0*se2 + ce0*ce2*se1) + w[2]*(ce0*se2 - ce2*se0*se1);
    double A12 = w[2]*ce0*ce1*ce2 - w[0]*ce2*se1 + w[1]*ce1*ce2*se0;
    double A13 = w[2]*(ce2*se0 - ce0*se1*se2) - w[1]*(ce0*ce2 + se0*se1*se2) - w[0]*ce1*se2;

    double A21 = -w[1]*(ce2*se0 - ce0*se1*se2) - w[2]*(ce0*ce2 + se0*se1*se2);
    double A22 = w[2]*ce0*ce1*se2 - w[0]*se1*se2 + w[1]*ce1*se0*se2;
    double A23 = w[2]*(se0*se2 + ce0*ce2*se1) - w[1]*(ce0*se2 - ce2*se0*se1) + w[0]*ce1*ce2;

    double A31 = w[1]*ce0*ce1 - w[2]*ce1*se0;
    double A32 = -w[0]*ce1 - w[2]*ce0*se1 -w[1]*se0*se1;
    double A33 = 0;
    
    Eigen::Matrix<double, 3, 3> A;
    A <<   A11,   A12,   A13,
           A21,   A22,   A23,
           A31,   A32,   A33;

    return A;
}

Eigen::Matrix<double, 3, 3> eulKinLinearization(const Eigen::Vector3d& e, const Eigen::Vector3d& w)
{
    double c_phi = cos(e[2]);
    double s_phi = sin(e[2]);
    double t_theta = tan(e[1]);
    double c_theta = cos(e[1]);
    double s_theta = sin(e[1]);

    double A11 = w[1]*c_phi*t_theta - w[2]*s_phi*t_theta;
    double A12 = w[2]*c_phi*(t_theta*t_theta + 1) + w[1]*s_phi*(t_theta*t_theta + 1);
    double A13 = 0;

    double A21 = - w[2]*c_phi - w[1]*s_phi;
    double A22 = 0;
    double A23 = 0;

    double A31 = (w[1]*c_phi)/c_theta - (w[2]*s_phi)/c_theta;
    double A32 = (w[2]*c_phi*s_theta)/c_theta*c_theta + (w[1]*s_phi*s_theta)/c_theta*c_theta;
    double A33 = 0;

    Eigen::Matrix<double, 3, 3> A;
    A <<   A11,   A12,   A13,
           A21,   A22,   A23,
           A31,   A32,   A33;

    return A;
}

// WGS

Eigen::Matrix<double, 3, 3> wgsToEnuRmat(double lat, double lon) // TODO not tested
{
    Eigen::Matrix<double, 3, 3> R_wgs_enu;
    R_wgs_enu <<
                -sin(lon)           ,  cos(lon)           ,  0       ,
                -sin(lat) * cos(lon), -sin(lat) * sin(lon),  cos(lat),
                 cos(lat) * cos(lon),  cos(lat) * sin(lon),  sin(lat);
                 
    return R_wgs_enu;
}

Eigen::Vector4d wgsToEnuQuat(double lat, double lon)
{
    double pi = 3.1415;
    Eigen::Vector4d qlat = quatFromEul(Eigen::Vector3d(pi/2 - lat, 0, 0));
    Eigen::Vector4d qlon = quatFromEul(Eigen::Vector3d(0, 0, pi/2 + lon));
                 
    Eigen::Vector4d res = quatMultiply(qlon, qlat);
    
    return res;
}
