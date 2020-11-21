#include "EKF4.h"

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <math.h>

#include <iostream>
#include <fstream>
#include <sstream>


EKF4::EKF4() :
	// consts
	O33(Eigen::Matrix<float, 3, 3>::Zero()),
	O34(Eigen::Matrix<float, 3, 4>::Zero()),
	O43(Eigen::Matrix<float, 4, 3>::Zero()),
	E33(Eigen::Matrix<float, 3, 3>::Identity())
{

}

void EKF4::initParams()
{
    // smoothed init
	_a_smoothed = Vector3(0.f, 0.f, 10.f);
	_w_smoothed = Vector3(0.f, 0.f, 0.f);
	_K_a_smoothed = 0.66f;
	_K_w_smoothed = 0.66f;
    
	// init Q
	Ekf4_state qDiag;
	qDiag << /*r*/ 1e-0f, 1e-0f, 1e-0f, /*v*/ 1e-1f, 1e-1f, 1e-1f, /*q*/ 1e-1f, 1e-1f, 1e-1f;
	_Q = qDiag.asDiagonal();

	// init P
	Ekf4_state pDiag = 50 * qDiag;
	_P = pDiag.asDiagonal();

	// init R_RV
	Vector6 rDiag_RV;
	rDiag_RV << /*r*/ 1e-3f, 1e-3f, 1e-3f, /*v*/ 1e-4f, 1e-4f, 1e-4f;
	_R_RV = rDiag_RV.asDiagonal();

	// init R_U
	Vector3 rDiag_U;
	rDiag_U <<  /*u*/ 1e-4f, 1e-4f, 1e-4f;
	_R_U = rDiag_U.asDiagonal();

	// init R_A
	Vector3 rDiag_A;
	rDiag_A << /*a*/ 1e-0f, 1e-0f, 1e-0f;
	_R_A = rDiag_A.asDiagonal();

	// init R_Q2
	Vector6 rDiag_Q2;
	rDiag_Q2 << 1e-6f, 1e-6f, 1e-6f, 1e-6f, 1e-6f, 1e-6f;
	_R_Q2 = rDiag_Q2.asDiagonal();


	//gpsAttachmentShift
	_drImuMaster << -0.28f, 0.0f, 0.1f;
	_drImuTarget << -0.28f, 0.0f, -0.48f;
	_drTargetMaster = _drImuMaster - _drImuTarget;
	_drSlave1 << 0.73f, 0.23f, 0.0f;
	_drSlave2 << 0.73f, -0.23f, 0.0f;
}

void EKF4::reset(const Vector3& r0)
{
    initParams();
    Vector3 v0(Vector3::Zero());
    Vector4 q0(1.0f, 0.0f, 0.0f, 0.0f);
    _X << r0, v0, q0;
}

void EKF4::calibSlavesWithSample(const Vector3& dr1, const Vector3& dr2)
{
    // assume slaves attached symmetric
    Vector3 d = dr1/2.0f + dr2/2.0f;
    Vector3 dh(d[0], d[1], 0.f);
    Vector4 q = quatBetweenVectors(dh, Vector3(1.f, 0.f, 0.f));
    _drSlave1 = quatRotate(q, dr1);
    _drSlave2 = quatRotate(q, dr2);
    //std::cout << _drSlave1 << std::endl;
    //std::cout << _drSlave2 << std::endl;
}

void EKF4::predict(float dt)
{
	/* state */
	Vector3 r0 = _X.segment(0, 3);
	Vector3 v0 = _X.segment(3, 3);
	Vector4 q0 = _X.segment(6, 4);

	Vector4 qw(0.0f, _w_smoothed[0], _w_smoothed[1], _w_smoothed[2]);
	Vector3 r1 = r0 + v0 * dt;
	Vector3 dq = quatRotate(q0, _w_smoothed).cross(v0) * dt;
	Vector3 v1 = v0 + quatRotate(q0, _w_smoothed).cross(v0) * dt;
	Vector4 q1 = q0 + 0.5 * quatMultiply(q0, qw) * dt;
	q1 = q1 / q1.norm();

	_X << r1, v1, q1;

	/* cov */
	Eigen::Matrix<float, 4, 4> Mqq_full = poissonEqLinearizationQ(_w_smoothed);
	Eigen::Matrix<float, 3, 3> Mqq = Mqq_full.block<3, 3>(1, 1);
	Eigen::Matrix<float, 4, 3> Mqw = poissonEqLinearizationW(q0, _w_smoothed);
	Eigen::Matrix<float, 3, 4> Mvq_full = -crossOperator(v0) * quatRotateLinearizationQ(q0, _w_smoothed);
	Eigen::Matrix<float, 3, 3> Mvq = Mvq_full.block<3, 3>(0, 1);

	Eigen::Matrix<float, 3, EKF4_STATE_DIM> Fr;
	Eigen::Matrix<float, 3, EKF4_STATE_DIM> Fv;
	Eigen::Matrix<float, 3, EKF4_STATE_DIM> Fq;

	Fr << O33, E33, O33;
	Fv << O33, O33, Mvq;
	Fq << O33, O33, Mqq;

	Eigen::Matrix<float, EKF4_STATE_DIM, EKF4_STATE_DIM> F;
	F << Fr, Fv, Fq;

	Eigen::Matrix<float, EKF4_STATE_DIM, EKF4_STATE_DIM> I(Eigen::Matrix<float, EKF4_STATE_DIM, EKF4_STATE_DIM>::Identity());
	Eigen::Matrix<float, EKF4_STATE_DIM, EKF4_STATE_DIM> PHI = I + F * dt;
	_P = PHI * _P * PHI + _Q;
}

void EKF4::correctRV(const Vector3& r, const Vector3& v)
{
    Vector6 rv;
    rv << r, v;
    correctRV(rv);
}

void EKF4::correctRV(const Vector6& rv)
{
	Vector3 r = _X.segment(0, 3);
	Vector3 v = _X.segment(3, 3);
	Vector4 q = _X.segment(6, 4);

	// mes model
	// Z[rgnns vgnns]
	Vector3 Zr = r + quatRotate(q, _drImuMaster);
	Vector3 Zv = v + quatRotate(q, _w_smoothed.cross(_drImuMaster));
	Vector6 Zx;
	Zx << Zr, Zv;
	Vector6 dz = rv - Zx;

	// H
	Eigen::Matrix<float, 3, 4> Zrq_full = quatRotateLinearizationQ(q, _drImuMaster);
	Eigen::Matrix<float, 3, 3> Zrq = Zrq_full.block<3, 3>(0, 1);
	Eigen::Matrix<float, 3, 4> Zvq_full = quatRotateLinearizationQ(q, _w_smoothed.cross(_drImuMaster));
	Eigen::Matrix<float, 3, 3> Zvq = Zvq_full.block<3, 3>(0, 1);

	Eigen::Matrix<float, 3, EKF4_STATE_DIM> H1;
	Eigen::Matrix<float, 3, EKF4_STATE_DIM> H2;
	Eigen::Matrix<float, 6, EKF4_STATE_DIM> H;
	H1 << E33, O33, Zrq;
	H2 << O33, E33, Zvq;
	H << H1, H2;


	// ordinary KF
	Eigen::Matrix<float, 6, 6> Rk = _R_RV + H * _P * H.transpose();
	Eigen::Matrix<float, EKF4_STATE_DIM, 6> K = _P * H.transpose() * (Rk.inverse());
	K = K;
	_P = _P - K * H * _P;
	Ekf4_state dx = K * dz;
	Ekf4_fullState dx_full;
	dx_full << dx.segment(0, 6), 0.f, dx.segment(6, 3);
	_X = _X + dx_full;
	_X.segment(6, 4).normalize();
}

void EKF4::correctQ2(const Vector3& dr1, const Vector3& dr2)
{
	Vector3 r = _X.segment(0, 3);
	Vector3 v = _X.segment(3, 3);
	Vector4 q = _X.segment(6, 4);

	// mes model
	Vector3 Z1 = quatRotate(q, _drSlave1);
	Vector3 Z2 = quatRotate(q, _drSlave2);
	Vector6 Zx;
	Zx << Z1, Z2;

	Vector6 Z;
	Z << dr1, dr2;
	Vector6 dz = Z - Zx;

	// H
	Eigen::Matrix<float, 3, 4> Zdr1_full = quatRotateLinearizationQ(q, _drSlave1);
	Eigen::Matrix<float, 3, 3> Zdr1 = Zdr1_full.block<3, 3>(0, 1);
	Eigen::Matrix<float, 3, 4> Zdr2_full = quatRotateLinearizationQ(q, _drSlave2);
	Eigen::Matrix<float, 3, 3> Zdr2 = Zdr2_full.block<3, 3>(0, 1);

	Eigen::Matrix<float, 3, EKF4_STATE_DIM> H1;
	Eigen::Matrix<float, 3, EKF4_STATE_DIM> H2;
	Eigen::Matrix<float, 6, EKF4_STATE_DIM> H;
	H1 << O33, O33, Zdr1;
	H2 << O33, O33, Zdr2;
	H << H1, H2;

	// ordinary KF
	Eigen::Matrix<float, 6, 6> Rk = _R_Q2 + H * _P * H.transpose();
	Eigen::Matrix<float, EKF4_STATE_DIM, 6> K = _P * H.transpose() * (Rk.inverse());
	_P = _P - K * H * _P;
	Ekf4_state dx = K * dz;
	Ekf4_fullState dx_full;
	dx_full << dx.segment(0, 6), 0.f, dx.segment(6, 3);
	_X = _X + dx_full;
	_X.segment(6, 4).normalize();
}


void EKF4::correctU(const Vector3& vMes)
{
	Vector3 r = _X.segment(0, 3);
	Vector3 v = _X.segment(3, 3);
	Vector4 q = _X.segment(6, 4);

	Vector3 ex(1.f, 0.f, 0.f);

	// mes model
	// Z[vgnns]
	Vector3 Zx = quatRotate(q, ex * v.norm()) + quatRotate(q, _w_smoothed.cross(_drTargetMaster));
	Vector3 dz = vMes - Zx;

	// H
	Eigen::Matrix<float, 3, 4> Zuq_full = quatRotateLinearizationQ(q, ex * v.norm() + _w_smoothed.cross(_drTargetMaster));
	Eigen::Matrix<float, 3, 3> Zuq = Zuq_full.block<3, 3>(0, 1);
	Eigen::Matrix<float, 3, 3> Zuv = Eigen::Matrix<float, 3, 3>::Zero();
	Zuv.row(0) = normVect3Linearization(v); // TODO  chek v is not zero
	Zuv = quatToMatrix(q) * Zuv;


	Eigen::Matrix<float, 3, EKF4_STATE_DIM> H;
	H << O33, Zuv, Zuq;

	// ordinary KF
	Eigen::Matrix<float, 3, 3> Rk = _R_U + H * _P * H.transpose();
	Eigen::Matrix<float, EKF4_STATE_DIM, 3> K = _P * H.transpose() * (Rk.inverse());
	_P = _P - K * H * _P;
	Ekf4_state dx = K * dz;
	Ekf4_fullState dx_full;
	dx_full << dx.segment(0, 6), 0.f, dx.segment(6, 3);
	_X = _X + dx_full;
	_X.segment(6, 4).normalize();
}

void EKF4::correctA(const Vector3& aMes)
{
	Vector3 r = _X.segment(0, 3);
	Vector3 v = _X.segment(3, 3);
	Vector4 q = _X.segment(6, 4);

	// TODO copy alg. from matlab.
}

Ekf4_fullState EKF4::getEstState()
{
	return _X;
}

Ekf4_fullState EKF4::getEstTargetState()
{
 	Vector3 r = _X.segment(0, 3);
	Vector3 v = _X.segment(3, 3);
	Vector4 q = _X.segment(6, 4);
 	Vector3 rTarget = r + quatRotate(q, _drImuTarget);
 	Vector3 vTarget = v + quatRotate(q, _drImuTarget.cross(_w_smoothed));
    Ekf4_fullState targetState;
    targetState << rTarget, vTarget, q;
    return targetState;
}

void EKF4::setImu(const Vector3& a, const Vector3& w)
{
	_a_smoothed = smooth(quatRotate(_qImuCalib, a), _a_smoothed, _K_a_smoothed);
	_w_smoothed = smooth(quatRotate(_qImuCalib, w), _w_smoothed, _K_w_smoothed);
}

void EKF4::setQImuCalib(const Vector4& q)
{
    _qImuCalib = q;
}

Vector3 EKF4::smooth(const Vector3& sample, const Vector3& smoothed, float K)
{
	Vector3 res = smoothed + K * (sample - smoothed);
	return res;
}

Vector3 EKF4::getDrSlave1()
{
    return _drSlave1;
}

Vector3 EKF4::getDrSlave2()
{
    return _drSlave2;
}
