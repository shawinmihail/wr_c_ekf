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
	qDiag << /*r*/ 1e-3f, 1e-3f, 1e-3f, /*v*/ 1e-4f, 1e-4f, 1e-4f, /*a*/ 1e-3f, 1e-3f, 1e-3f, /*q*/ 1e-8f, 1e-8f, 1e-8f, 1e-8f, /*w*/ 1e-8f, 1e-8f, 1e-8f;
	_Q = qDiag.asDiagonal();

	// init P
	EkfStateVector pDiag = 25 * qDiag;
	Eigen::Matrix<float, SREKF_STATE_DIM, SREKF_STATE_DIM> P = pDiag.asDiagonal();
	Eigen::LLT<Eigen::Matrix<float, SREKF_STATE_DIM, SREKF_STATE_DIM>> cdP(P);
	_sqrtP = cdP.matrixL();

	// init R_pv
	Vector6 rDiag_pv;
	rDiag_pv << /*r*/ 1e-1f, 1e-1f, 1e-1f, /*v*/ 1e-2f, 1e-2f, 1e-2f;
	Eigen::Matrix<float, 6, 6> R_pv = rDiag_pv.asDiagonal();
	Eigen::LLT<Eigen::Matrix<float, 6, 6>> cdR_pv(R_pv);
	_sqrtR_pv = cdR_pv.matrixL();

	// init R_v
	Vector3 rDiag_v;
	rDiag_v <<  /*v*/ 1e-2f, 1e-2f, 1e-2f;
	Eigen::Matrix<float, 3, 3> R_v = rDiag_v.asDiagonal();
	Eigen::LLT<Eigen::Matrix<float, 3, 3>> cdR_v(R_v);
	_sqrtR_v = cdR_v.matrixL();

	// init R_a
	Vector3 rDiag_a;
	rDiag_a << /*a*/ 1e-2f, 1e-2f, 1e-2f;
	Eigen::Matrix<float, 3, 3> R_a = rDiag_a.asDiagonal();
	Eigen::LLT<Eigen::Matrix<float, 3, 3>> cdR_a(R_a);
	_sqrtR_a = cdR_a.matrixL();

	//gpsAttachmentShift
	_gpsAttachmentShift << 0.3f, -0.3f, 0.9f;

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
	triaArg << Enn + F * dt * _sqrtP, _Q * sqrt(dt);

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

	//%% ordinary K, H
	//	% P = sqrtP * sqrtP';
	//	% R = sqrtR * sqrtR';
	//	% Rk = R + H * P * H';
	//	% K = P * H' * (Rk)^-1;
	//	% P = P - K * H * P;
	//% sqrtP = chol(P, 'lower');
	//% X1 = X + K * dz;

	// square - root K, H
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

EkfStateVector SREKF::getEstState()
{
	return _X;
}

