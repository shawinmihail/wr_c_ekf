#include "EKF4_WGS.h"

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <math.h>

#include <iostream>
#include <fstream>
#include <sstream>


EKF4_WGS::EKF4_WGS() :
	// consts
	O33(Eigen::Matrix<double, 3, 3>::Zero()),
	O34(Eigen::Matrix<double, 3, 4>::Zero()),
	O43(Eigen::Matrix<double, 4, 3>::Zero()),
	E33(Eigen::Matrix<double, 3, 3>::Identity()),
	_daImuCalib(Eigen::Vector3d::Zero()),
	_dwImuCalib(Eigen::Vector3d::Zero())
{

}

/* -> INIT */

void EKF4_WGS::initParams()
{
    // smoothed init
	_a_smoothed = Eigen::Vector3d(0, 0, 10);
	_w_smoothed = Eigen::Vector3d(0, 0, 0);
	_K_a_smoothed = 0.33;
	_K_w_smoothed = 0.33;
    
	// init Q
	Eigen::Matrix<double, 10, 1> qDiag;
	qDiag << /*r*/ 1e-0, 1e-0, 1e-0, /*v*/ 1e-1, 1e-1, 1e-1, /*q*/ 1e-4, 1e-4, 1e-4, 1e-4;
	_Q = qDiag.asDiagonal();

	// init P
	Eigen::Matrix<double, 10, 1> pDiag = 50 * qDiag;
	_P = pDiag.asDiagonal();

	// init R_RV
	Eigen::Matrix<double, 6, 1> rDiag_RV;
	rDiag_RV << /*r*/ 1e-3, 1e-3, 1e-3, /*v*/ 1e-4, 1e-4, 1e-4;
	_R_RV = rDiag_RV.asDiagonal();

	// init R_Q2
	Eigen::Matrix<double, 6, 1> rDiag_Q2;
	rDiag_Q2 << 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4;
	_R_Q2 = rDiag_Q2.asDiagonal();
}

bool EKF4_WGS::initState(const Eigen::Vector3d& r0, const Eigen::Vector3d& s1, const Eigen::Vector3d& s2, double slaves_accuracy)
{
    Eigen::Vector3d v0(Eigen::Vector3d::Zero());
    Eigen::Vector4d q0(1.0, 0.0, 0.0, 0.0);
    bool is_ok = attitudeWithBaseLineNumeric(s1, s2, slaves_accuracy, q0);
    
    if (is_ok)
    {
        _X << r0, v0, q0;
    }
    
    return is_ok;
}

/* INIT <- */

/* -> EKF */

void EKF4_WGS::predict(double dt)
{
	/* state */
	Eigen::Vector3d r0 = _X.segment(0, 3);
	Eigen::Vector3d v0 = _X.segment(3, 3);
	Eigen::Vector4d q0 = _X.segment(6, 4);

	Eigen::Vector4d qw(0.0, _w_smoothed[0], _w_smoothed[1], _w_smoothed[2]);
	Eigen::Vector3d r1 = r0 + v0 * dt;
	Eigen::Vector3d v1 = v0 + quatRotate(q0, _w_smoothed).cross(v0) * dt;
	Eigen::Vector4d q1 = q0 + 0.5 * quatMultiply(q0, qw) * dt;
	q1 = q1 / q1.norm();

	_X << r1, v1, q1;

	/* cov */
	Eigen::Matrix<double, 4, 4> Mqq_full = poissonEqLinearizationQ(_w_smoothed);
	Eigen::Matrix<double, 4, 3> Mqw = poissonEqLinearizationW(q0, _w_smoothed);
	Eigen::Matrix<double, 3, 4> Mvq_full = -crossOperator(v0) * quatRotateLinearizationQ(q0, _w_smoothed);

	Eigen::Matrix<double, 3, 10> Fr;
	Eigen::Matrix<double, 3, 10> Fv;
	Eigen::Matrix<double, 4, 10> Fq;

	Fr << O33, E33, O34;
	Fv << O33, O33, Mvq_full;
	Fq << O43, O43, Mqq_full;

	Eigen::Matrix<double, 10, 10> F;
	F << Fr, Fv, Fq;

	Eigen::Matrix<double, 10, 10> I(Eigen::Matrix<double, 10, 10>::Identity());
	Eigen::Matrix<double, 10, 10> PHI = I + F * dt;
	_P = PHI * _P * PHI.transpose() + _Q;
}

void EKF4_WGS::correctRV(const Eigen::Vector3d& r_mes, const Eigen::Vector3d& v_mes)
{

    Eigen::Vector3d r = _X.segment(0, 3);
	Eigen::Vector3d v = _X.segment(3, 3);
	Eigen::Vector4d q = _X.segment(6, 4);
    Eigen::Matrix<double, 6, 1> rv_mes;
    rv_mes << r_mes, v_mes;

	// mes model
	// Z[rgnns vgnns]
	Eigen::Vector3d Zr = r + quatRotate(q, _drImuMaster);
	Eigen::Vector3d Zv = v + quatRotate(q, _w_smoothed.cross(_drImuMaster));
	Eigen::Matrix<double, 6, 1> Zx;
	Zx << Zr, Zv;
	Eigen::Matrix<double, 6, 1> dz = rv_mes - Zx;

	// H
	Eigen::Matrix<double, 3, 4> Zrq_full = quatRotateLinearizationQ(q, _drImuMaster);
	Eigen::Matrix<double, 3, 4> Zvq_full = quatRotateLinearizationQ(q, _w_smoothed.cross(_drImuMaster));

	Eigen::Matrix<double, 3, 10> H1;
	Eigen::Matrix<double, 3, 10> H2;
	Eigen::Matrix<double, 6, 10> H;
	H1 << E33, O33, Zrq_full;
	H2 << O33, E33, Zvq_full;
	H << H1, H2;


	// ordinary KF
	Eigen::Matrix<double, 6, 6> Rk = _R_RV + H * _P * H.transpose();
	Eigen::Matrix<double, 10, 6> K = _P * H.transpose() * (Rk.inverse());
	_P = _P - K * H * _P;
	Eigen::Matrix<double, 10, 1> dx = K * dz;
	_X = _X + dx;
	_X.segment(6, 4).normalize();
}

void EKF4_WGS::correctQ2(const Eigen::Vector3d& dr1, const Eigen::Vector3d& dr2)
{
	Eigen::Vector3d r = _X.segment(0, 3);
	Eigen::Vector3d v = _X.segment(3, 3);
	Eigen::Vector4d q = _X.segment(6, 4);

	// mes model
	Eigen::Vector3d Z1 = quatRotate(q, _drSlave1);
	Eigen::Vector3d Z2 = quatRotate(q, _drSlave2);
	Eigen::Matrix<double, 6, 1> Zx;
	Zx << Z1, Z2;

	Eigen::Matrix<double, 6, 1> Z;
	Z << dr1, dr2;
	Eigen::Matrix<double, 6, 1> dz = Z - Zx;

	// H
	Eigen::Matrix<double, 3, 4> Zdr1_full = quatRotateLinearizationQ(q, _drSlave1);
	Eigen::Matrix<double, 3, 4> Zdr2_full = quatRotateLinearizationQ(q, _drSlave2);

	Eigen::Matrix<double, 3, 10> H1;
	Eigen::Matrix<double, 3, 10> H2;
	Eigen::Matrix<double, 6, 10> H;
	H1 << O33, O33, Zdr1_full;
	H2 << O33, O33, Zdr2_full;
	H << H1, H2;

	// ordinary KF
	Eigen::Matrix<double, 6, 6> Rk = _R_Q2 + H * _P * H.transpose();
	Eigen::Matrix<double, 10, 6> K = _P * H.transpose() * (Rk.inverse());
	_P = _P - K * H * _P;
	Eigen::Matrix<double, 10, 1> dx = K * dz;
	_X = _X + dx;
	_X.segment(6, 4).normalize();
}

/* INIT <- */

/* -> */

bool EKF4_WGS::attitudeWithBaseLineNumeric(const Eigen::Vector3d& s1, const Eigen::Vector3d& s2, double slaves_accuracy, Eigen::Vector4d& q_out)
{
    bool is_ok = false;
    
    double eps = 0.5;
    for (int i = 0; i < 99; i++)
    {
        Eigen::Matrix<double, 6, 1> d;
        d << quatRotate(q_out, _drSlave1) - s1,  quatRotate(q_out, _drSlave2) - s2;
        Eigen::Matrix<double, 3, 4> Jq1 = quatRotateLinearizationQ(q_out, _drSlave1);
        Eigen::Matrix<double, 3, 4> Jq2 = quatRotateLinearizationQ(q_out, _drSlave2);
        Eigen::MatrixXd Jq(6, 4);
        Jq << Jq1, Jq2;
        Eigen::MatrixXd A = Jq.completeOrthogonalDecomposition().pseudoInverse();
        
        q_out = q_out - eps * A * d;
        q_out = q_out / q_out.norm();
        double y = d.norm();
                
        if (y < slaves_accuracy*2)
        {
            is_ok = true;
            break;
        }
    }
    
    return is_ok;
}


Eigen::Vector3d EKF4_WGS::smooth(const Eigen::Vector3d& sample, const Eigen::Vector3d& smoothed, double K)
{
	Eigen::Vector3d res = smoothed + K * (sample - smoothed);
	return res;
}

/* <- */

/* -> SET */

void EKF4_WGS::setSlavesCalib(const Eigen::Vector3d& slave1, const Eigen::Vector3d& slave2)
{
    _drSlave1 = slave1;
    _drSlave2 = slave2;
}

void EKF4_WGS::setSensorsGeometry(const Eigen::Vector3d& drImuMaster, const Eigen::Vector3d& drImuTarget)
{
    _drImuMaster = drImuMaster;
    _drImuTarget = drImuTarget;
    _drTargetMaster = _drImuMaster - _drImuTarget;
};

void EKF4_WGS::setImu(const Eigen::Vector3d& a, const Eigen::Vector3d& w)
{
	_a_smoothed = smooth(quatRotate(_qImuCalib, a - _daImuCalib), _a_smoothed, _K_a_smoothed);
	_w_smoothed = smooth(quatRotate(_qImuCalib, w - _dwImuCalib), _w_smoothed, _K_w_smoothed);
}

void EKF4_WGS::setQImuCalib(const Eigen::Vector4d& q)
{
    _qImuCalib = q;
}

void EKF4_WGS::setBiasesImuCalib(const Eigen::Vector3d& da, const Eigen::Vector3d& dw)
{
    _daImuCalib = da;
    _dwImuCalib = dw;
}

/* SET <- */

/* -> GET */

Eigen::Matrix<double, 10, 1> EKF4_WGS::getEstState()
{
	return _X;
}

Eigen::Matrix<double, 10, 1> EKF4_WGS::getEstTargetState()
{
 	Eigen::Vector3d r = _X.segment(0, 3);
	Eigen::Vector3d v = _X.segment(3, 3);
	Eigen::Vector4d q = _X.segment(6, 4);
 	Eigen::Vector3d rTarget = r + quatRotate(q, _drImuTarget);
 	Eigen::Vector3d vTarget = v + quatRotate(q, _drImuTarget.cross(_w_smoothed));
    Eigen::Matrix<double, 10, 1> targetState;
    targetState << rTarget, vTarget, q;
    return targetState;
}

Eigen::Vector3d EKF4_WGS::getDrSlave1()
{
    return _drSlave1;
}

Eigen::Vector3d EKF4_WGS::getDrSlave2()
{
    return _drSlave2;
}

/* GET <- */
