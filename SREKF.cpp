#include "SREKF.h"

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <math.h>

#include <iostream>
#include <fstream>
#include <sstream>

template <int a, int b>
Eigen::Matrix<float, a, a> getLowTriang(const Eigen::Matrix<float, a, b>& A)
{
	Eigen::Matrix<float, b, a> A_transposed = A.transpose();
	Eigen::HouseholderQR<Eigen::Matrix<float, b, a>> qrHolder;
	qrHolder.compute(A_transposed);
	Eigen::Matrix<float, b, a> R = qrHolder.matrixQR().triangularView<Eigen::Upper>();
	Eigen::Matrix<float, a, b> L = R.transpose();
	Eigen::Matrix<float, a, a> cutedL = L.block(0, 0, a, a);

	return -cutedL;
}

template <int a>
Eigen::Matrix<float, a, a> cholUpdate(const Eigen::Matrix<float, a, a>& A)
{
	Eigen::LLT<Eigen::Matrix<float, a, a>> cdA(A);
	return cdA.matrixL();
}

SREKF::SREKF() :
	// consts
	O33(Eigen::Matrix<float, 3, 3>::Zero()),
	O34(Eigen::Matrix<float, 3, 4>::Zero()),
	O43(Eigen::Matrix<float, 4, 3>::Zero()),
	E33(Eigen::Matrix<float, 3, 3>::Identity())
{
	_X = EkfStateVector::Zero();
	_X(9) = 1.f;

	// init Q
	EkfStateVector qDiag;
	qDiag << /*r*/ 1e-2f, 1e-2f, 1e-2f, /*v*/ 1e-2f, 1e-2f, 1e-2f, /*a*/ 1e-2f, 1e-2f, 1e-2f, /*q*/ 1e-4f, 1e-4f, 1e-4f, 1e-4f, /*w*/ 1e-2f, 1e-2f, 1e-2f;
	_Q = qDiag.asDiagonal();

	// init P
	EkfStateVector pDiag = 10 * qDiag;
	Eigen::Matrix<float, SREKF_STATE_DIM, SREKF_STATE_DIM> P = pDiag.asDiagonal();
	_sqrtP = cholUpdate<SREKF_STATE_DIM>(P);

	// init R_pv
	Vector6 rDiag_pv;
	rDiag_pv << /*r*/ 0.144f, 0.144f, 0.144f, /*v*/ 0.0025f, 0.0025f, 0.0025f;
	Eigen::Matrix<float, 6, 6> R_pv = rDiag_pv.asDiagonal();
	_sqrtR_pv = cholUpdate<6>(R_pv);

	// init R_v
	Vector3 rDiag_v;
	rDiag_v <<  /*v*/ 0.0025f, 0.0025f, 0.0025f;
	Eigen::Matrix<float, 3, 3> R_v = rDiag_v.asDiagonal();
	_sqrtR_v = cholUpdate<3>(R_v);

	// init R_a
	Vector3 rDiag_a;
	rDiag_a << /*a*/ 0.0025f, 0.0025f, 0.0025f;
	Eigen::Matrix<float, 3, 3> R_a = rDiag_a.asDiagonal();
	_sqrtR_a = cholUpdate<3>(R_a);

	// init R_p3
	Vector6 rDiag_p3;
	rDiag_p3 <<  0.144f, 0.144f, 0.144f, 0.144f, 0.144f, 0.144f;
	Eigen::Matrix<float, 6, 6> R_p3 = rDiag_p3.asDiagonal();
	_sqrtR_p3 = cholUpdate<6>(R_p3);

	// init R_z
	Eigen::Matrix<float, 12, 1> rDiag_z;
	rDiag_z << /*r*/ 15e-3f, 15e-3f, 15e-3f, /*v*/ 3e-3f, 3e-3f, 3e-3f, /*v*/ 3e-3f, 3e-3f, 3e-3f, /*a*/ 3e-3f, 3e-3f, 3e-3f;
	Eigen::Matrix<float, 12, 12> R_z = rDiag_z.asDiagonal();
	_sqrtR_z = cholUpdate<12>(R_z);

	//gpsAttachmentShift
	_gpsAttachmentShift << 0.3f, -0.3f, 0.9f;
	_gpsSlave1 << 0.3f, 0.2f, 0.1f;
	_gpsSlave2 << 0.3f, -0.2f, 0.1f;
}

void SREKF::predictImu(const Vector3& aMes, const Vector3& wMes, float dt)
{

	/* state */
	Vector3 g;
	g << 0.f, 0.f, -10.f;
	Vector3 r0 = _X.segment(0, 3);
	Vector3 v0 = _X.segment(3, 3);
	Vector4 q0 = _X.segment(9, 4);
	Vector3 w0 = wMes;
	Vector3 a0 = quatRotate(q0, aMes) + g;

	Vector4 qw0(0.0f, w0[0], w0[1], w0[2]);
	Vector3 r1 = r0 + v0 * dt;
	Vector3 v1 = v0 + a0 * dt;
	Vector3 a1 = a0;
	Vector4 q1 = q0 + 0.5 * quatMultiply(q0, qw0) * dt;
	//Vector3 q1 = q_next / norm(q_next);
	Vector3 w1 = w0;

	_X << r1, v1, a1, q1, w1;

	/* cov */
	Eigen::Matrix<float, 3, 4> Maq = quatRotateLinearizationQ(q0, aMes);
	Eigen::Matrix<float, 4, 4> Mqq = poissonEqLinearizationQ(wMes);
	Eigen::Matrix<float, 4, 3> Mqw = poissonEqLinearizationW(q0, wMes);

	Eigen::Matrix<float, 3, SREKF_STATE_DIM> Fr;
	Eigen::Matrix<float, 3, SREKF_STATE_DIM> Fv;
	Eigen::Matrix<float, 3, SREKF_STATE_DIM> Fa;
	Eigen::Matrix<float, 4, SREKF_STATE_DIM> Fq;
	Eigen::Matrix<float, 3, SREKF_STATE_DIM> Fw;

	Fr << O33, E33, O33, O34, O33;
	Fv << O33, O33, E33, Maq, O33;
	Fa << O33, O33, O33, O34, O33;
	Fq << O43, O43, O43, Mqq, Mqw;
	Fw << O33, O33, O33, O34, O33;

	Eigen::Matrix<float, SREKF_STATE_DIM, SREKF_STATE_DIM> F;
	F << Fr, Fv, Fa, Fq, Fw;

	Eigen::Matrix<float, SREKF_STATE_DIM, SREKF_STATE_DIM> Enn(Eigen::Matrix<float, SREKF_STATE_DIM, SREKF_STATE_DIM>::Identity());

	Eigen::Matrix<float, SREKF_STATE_DIM, 2 * SREKF_STATE_DIM> triaArg;
	triaArg << (Enn + F * dt) * _sqrtP, _Q * sqrt(dt);

	_sqrtP = getLowTriang<SREKF_STATE_DIM, 2 * SREKF_STATE_DIM>(triaArg);
}

void SREKF::correctPv(const Vector6& pv)
{
	Vector3 r = _X.segment(0, 3);
	Vector3 v = _X.segment(3, 3);
	Vector3 a = _X.segment(6, 3);
	Vector4 q = _X.segment(9, 4);
	Vector3 w = _X.segment(13, 3);

	// mes model
	// Z[rgnns vgnns]
	Vector3 Zr = r + quatRotate(q, _gpsAttachmentShift);
	Vector3 Zv = v + quatRotate(q, w.cross(_gpsAttachmentShift));
	Vector6 Zx;
	Zx << Zr, Zv;
	Vector6 dz = pv - Zx;

	// H
	Eigen::Matrix<float, 3, 4> Zrq = quatRotateLinearizationQ(q, _gpsAttachmentShift);
	Eigen::Matrix<float, 3, 4> Zvq = quatRotateLinearizationQ(q, w.cross(_gpsAttachmentShift));
	Eigen::Matrix<float, 3, 3> Zvw = quatToMatrix(q) * crossOperator(-_gpsAttachmentShift);

	Eigen::Matrix<float, 3, SREKF_STATE_DIM> H1;
	Eigen::Matrix<float, 3, SREKF_STATE_DIM> H2;
	Eigen::Matrix<float, 6, SREKF_STATE_DIM> H;
	H1 << E33, O33, O33, Zrq, O33;
	H2 << O33, E33, O33, Zvq, Zvw;
	H << H1, H2;

	// square-root
	Eigen::Matrix<float, SREKF_STATE_DIM + 6, SREKF_STATE_DIM + 6> triaArg;
	Eigen::Matrix<float, 6, SREKF_STATE_DIM + 6> triaArg1;
	Eigen::Matrix<float, SREKF_STATE_DIM, SREKF_STATE_DIM + 6> triaArg2;
	triaArg1 << _sqrtR_pv, H * _sqrtP;

	triaArg2 << Eigen::Matrix<float, SREKF_STATE_DIM, 6>::Zero(), _sqrtP;
	triaArg << triaArg1, triaArg2;

	Eigen::Matrix<float, SREKF_STATE_DIM + 6, SREKF_STATE_DIM + 6> M = getLowTriang<SREKF_STATE_DIM + 6, SREKF_STATE_DIM + 6>(triaArg);
	Eigen::Matrix<float, 6, 6> sqrtRk = M.block<6, 6>(0, 0);
	Eigen::Matrix<float, SREKF_STATE_DIM, 6> K = M.block<SREKF_STATE_DIM, 6>(6, 0);
	_sqrtP = M.block<SREKF_STATE_DIM, SREKF_STATE_DIM>(6, 6);

	_X = _X + K * (sqrtRk.transpose().inverse())*dz;
	_X.segment(9, 4).normalize();
}

void SREKF::correctP3(const Vector3& dr1, const Vector3& dr2)
{
	Vector3 r = _X.segment(0, 3);
	Vector3 v = _X.segment(3, 3);
	Vector3 a = _X.segment(6, 3);
	Vector4 q = _X.segment(9, 4);
	Vector3 w = _X.segment(13, 3);

	// mes model
	Vector3 Z1 = quatRotate(q, _gpsSlave1);
	Vector3 Z2 = quatRotate(q, _gpsSlave2);
	Vector6 Zx;
	Zx << Z1, Z2;

	Vector6 Z;
	Z << dr1, dr2;
	Vector6 dz = Z - Zx;

	// H
	Eigen::Matrix<float, 3, 4> Zdr1 = quatRotateLinearizationQ(q, _gpsSlave1);
	Eigen::Matrix<float, 3, 4> Zdr2 = quatRotateLinearizationQ(q, _gpsSlave2);

	Eigen::Matrix<float, 3, SREKF_STATE_DIM> H1;
	Eigen::Matrix<float, 3, SREKF_STATE_DIM> H2;
	Eigen::Matrix<float, 6, SREKF_STATE_DIM> H;
	H1 << O33, O33, O33, Zdr1, O33;
	H2 << O33, O33, O33, Zdr2, O33;
	H << H1, H2;

	// square-root
	Eigen::Matrix<float, SREKF_STATE_DIM + 6, SREKF_STATE_DIM + 6> triaArg;
	Eigen::Matrix<float, 6, SREKF_STATE_DIM + 6> triaArg1;
	Eigen::Matrix<float, SREKF_STATE_DIM, SREKF_STATE_DIM + 6> triaArg2;
	triaArg1 << _sqrtR_p3, H * _sqrtP;

	triaArg2 << Eigen::Matrix<float, SREKF_STATE_DIM, 6>::Zero(), _sqrtP;
	triaArg << triaArg1, triaArg2;

	Eigen::Matrix<float, SREKF_STATE_DIM + 6, SREKF_STATE_DIM + 6> M = getLowTriang<SREKF_STATE_DIM + 6, SREKF_STATE_DIM + 6>(triaArg);
	Eigen::Matrix<float, 6, 6> sqrtRk = M.block<6, 6>(0, 0);
	Eigen::Matrix<float, SREKF_STATE_DIM, 6> K = M.block<SREKF_STATE_DIM, 6>(6, 0);
	_sqrtP = M.block<SREKF_STATE_DIM, SREKF_STATE_DIM>(6, 6);

	_X = _X + K * (sqrtRk.transpose().inverse())*dz;
	_X.segment(9, 4).normalize();
}

void SREKF::correctZ(const Vector3& pMes, const Vector3& vMes, const Vector3& aMes)
{
	Eigen::Matrix<float, 12, 1> Z;
	Z << pMes, vMes, vMes, aMes;

	// mes model
	Vector3 r = _X.segment(0, 3);
	Vector3 v = _X.segment(3, 3);
	Vector3 a = _X.segment(6, 3);
	Vector4 q = _X.segment(9, 4);
	Vector3 w = _X.segment(13, 3);

	Vector3 ex(1.f, 0.f, 0.f);
	Vector3 g(0.f, 0.f, -10.f);

	Vector3 Zr = r + quatRotate(q, _gpsAttachmentShift);
	Vector3 Zv = v + quatRotate(q, w.cross(_gpsAttachmentShift));
	Vector3 Zu = quatRotate(q, ex * v.norm()) + quatRotate(q, w.cross(_gpsAttachmentShift));
	Vector3 Za = quatRotate(quatInverse(q), a - g);

	Eigen::Matrix<float, 12, 1> Zx;
	Zx << Zr, Zu, Zv, Za;
	Eigen::Matrix<float, 12, 1> dz = Z - Zx;

	// H
	Eigen::Matrix<float, 3, 4> Zrq = quatRotateLinearizationQ(q, _gpsAttachmentShift);
	Eigen::Matrix<float, 3, 4> Zvq = quatRotateLinearizationQ(q, w.cross(_gpsAttachmentShift));
	Eigen::Matrix<float, 3, 3> Zvw = quatToMatrix(q) * crossOperator(-_gpsAttachmentShift);

	Eigen::Matrix<float, 3, 4> Zuq = quatRotateLinearizationQ(q, ex * v.norm()) + Zvq;
	Eigen::Matrix<float, 3, 3> Zuv = Eigen::Matrix<float, 3, 3>::Zero();
	Zuv.row(0) = normVect3Linearization(v); // !!! v not zero chek add
	Zuv = quatToMatrix(q) * Zuv;
	Eigen::Matrix<float, 3, 3> Zuw = Zvw;

	Eigen::Matrix<float, 3, 4> Zaq = quatRotateLinearizationQ(quatInverse(q), a - g);
	Zaq.block<3, 3>(0, 1) = -Zaq.block<3, 3>(0, 1);
	Eigen::Matrix<float, 3, 3> Zaa = quatToMatrix(quatInverse(q));

	Eigen::Matrix<float, 3, SREKF_STATE_DIM> Hr;
	Eigen::Matrix<float, 3, SREKF_STATE_DIM> Hv;
	Eigen::Matrix<float, 3, SREKF_STATE_DIM> Hu;
	Eigen::Matrix<float, 3, SREKF_STATE_DIM> Ha;
	Eigen::Matrix<float, 12, SREKF_STATE_DIM> H;
	Hr << E33, O33, O33, Zrq, O33;
	Hv << O33, E33, O33, Zvq, Zvw;
	Hu << O33, Zuv, O33, Zuq, Zuw;
	Ha << O33, O33, Zaa, Zaq, O33;
	H << Hr, Hu, Hv, Ha;

	// ordinary
	Eigen::Matrix<float, SREKF_STATE_DIM, SREKF_STATE_DIM> P = _sqrtP * _sqrtP.transpose();
	Eigen::Matrix<float, 12, 12> R = _sqrtR_z * _sqrtR_z.transpose();
	Eigen::Matrix<float, 12, 12> Rk = R + H * P * H.transpose();
	Eigen::Matrix<float, SREKF_STATE_DIM, 12> K = P * H.transpose() * (Rk.inverse());
	P = P - K * H * P;
	_sqrtP = cholUpdate<SREKF_STATE_DIM>(P);
	_X = _X + K * dz;
	_X.segment(9, 4).normalize();
}

void SREKF::correctV(const Vector3& vMes)
{
	Vector3 r = _X.segment(0, 3);
	Vector3 v = _X.segment(3, 3);
	Vector3 a = _X.segment(6, 3);
	Vector4 q = _X.segment(9, 4);
	Vector3 w = _X.segment(13, 3);

	Vector3 ex(1.f, 0.f, 0.f);

	// mes model
	// Z[vgnns]
	Vector3 Zx = quatRotate(q, ex * v.norm()) + quatRotate(q, w.cross(_gpsAttachmentShift));
	Vector3 dz = vMes - Zx;

	// H
	Eigen::Matrix<float, 3, 4> Zuq = quatRotateLinearizationQ(q, ex * v.norm() + w.cross(_gpsAttachmentShift));
	Eigen::Matrix<float, 3, 3> Zuv = Eigen::Matrix<float, 3, 3>::Zero();
	Zuv.row(0) = normVect3Linearization(v); // !!! v not zero chek add
	Zuv = quatToMatrix(q) * Zuv;
	Eigen::Matrix<float, 3, 3> Zuw = quatToMatrix(q) * crossOperator(-_gpsAttachmentShift);

	Eigen::Matrix<float, 3, SREKF_STATE_DIM> H;
	H << O33, Zuv, O33, Zuq, Zuw;

	// square-root
	Eigen::Matrix<float, SREKF_STATE_DIM + 3, SREKF_STATE_DIM + 3> triaArg;
	Eigen::Matrix<float, 3, SREKF_STATE_DIM + 3> triaArg1;
	Eigen::Matrix<float, SREKF_STATE_DIM, SREKF_STATE_DIM + 3> triaArg2;
	triaArg1 << _sqrtR_v, H * _sqrtP;

	triaArg2 << Eigen::Matrix<float, SREKF_STATE_DIM, 3>::Zero(), _sqrtP;
	triaArg << triaArg1, triaArg2;

	Eigen::Matrix<float, SREKF_STATE_DIM + 3, SREKF_STATE_DIM + 3> M = getLowTriang<SREKF_STATE_DIM + 3, SREKF_STATE_DIM + 3>(triaArg);
	Eigen::Matrix<float, 3, 3> sqrtRk = M.block<3, 3>(0, 0);
	Eigen::Matrix<float, SREKF_STATE_DIM, 3> K = M.block<SREKF_STATE_DIM, 3>(3, 0);
	_sqrtP = M.block<SREKF_STATE_DIM, SREKF_STATE_DIM>(3, 3);

	_X = _X + K * (sqrtRk.transpose().inverse()) * dz;
	_X.segment(9, 4).normalize();
}

void SREKF::correctA(const Vector3& aMes)
{
	Vector3 g(0.f, 0.f, -10.f);

	Vector3 r = _X.segment(0, 3);
	Vector3 v = _X.segment(3, 3);
	Vector3 a = _X.segment(6, 3);
	Vector4 q = _X.segment(9, 4);
	Vector3 w = _X.segment(13, 3);

	Vector3 ex(1.f, 0.f, 0.f);

	// mes model
	// Z[vgnns]
	Vector3 Zx = quatRotate(quatInverse(q), a - g);
	Vector3 dz = aMes - Zx;

	// H
	Eigen::Matrix<float, 3, 4> Zaq = quatRotateLinearizationQ(quatInverse(q), a - g);
	Zaq.block<3, 3>(0, 1) = -Zaq.block<3, 3>(0, 1);
	Eigen::Matrix<float, 3, 3> Zaa = quatToMatrix(quatInverse(q));

	Eigen::Matrix<float, 3, SREKF_STATE_DIM> H;
	H << O33, O33, Zaa, Zaq, O33;

	// square-root
	Eigen::Matrix<float, SREKF_STATE_DIM + 3, SREKF_STATE_DIM + 3> triaArg;
	Eigen::Matrix<float, 3, SREKF_STATE_DIM + 3> triaArg1;
	Eigen::Matrix<float, SREKF_STATE_DIM, SREKF_STATE_DIM + 3> triaArg2;
	triaArg1 << _sqrtR_a, H* _sqrtP;

	triaArg2 << Eigen::Matrix<float, SREKF_STATE_DIM, 3>::Zero(), _sqrtP;
	triaArg << triaArg1, triaArg2;

	Eigen::Matrix<float, SREKF_STATE_DIM + 3, SREKF_STATE_DIM + 3> M = getLowTriang<SREKF_STATE_DIM + 3, SREKF_STATE_DIM + 3>(triaArg);
	Eigen::Matrix<float, 3, 3> sqrtRk = M.block<3, 3>(0, 0);
	Eigen::Matrix<float, SREKF_STATE_DIM, 3> K = M.block<SREKF_STATE_DIM, 3>(3, 0);
	_sqrtP = M.block<SREKF_STATE_DIM, SREKF_STATE_DIM>(3, 3);

	_X = _X + K * (sqrtRk.transpose().inverse()) * dz;
	_X.segment(9, 4).normalize();
}

EkfStateVector SREKF::getEstState()
{
	return _X;
}

