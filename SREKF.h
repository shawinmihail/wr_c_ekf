#pragma once
#include "Definitions.h"
#include "Utils.h"

// state X = [r v a q w]
#define SREKF_STATE_DIM 16
typedef Eigen::Matrix<float, SREKF_STATE_DIM, 1> EkfStateVector;


class SREKF
{
public:
	SREKF();
	void predictImu(const Vector3& aMes, const Vector3& wMes, float dt);
	void SREKF::correctZ(const Vector3& p, const Vector3& v, const Vector3& a);
	void correctPv(const Vector6& pv);
	void correctV(const Vector3& v);
	void correctA(const Vector3& a);
	EkfStateVector getEstState();

private:
	EkfStateVector _X;
	Eigen::Matrix<float, SREKF_STATE_DIM, SREKF_STATE_DIM> _Q;
	Eigen::Matrix<float, SREKF_STATE_DIM, SREKF_STATE_DIM> _sqrtP;

	Eigen::Matrix<float, 6, 6> _sqrtR_pv;
	Eigen::Matrix<float, 3, 3> _sqrtR_v;
	Eigen::Matrix<float, 3, 3> _sqrtR_a;
	Eigen::Matrix<float, 12, 12> _sqrtR_z;

	Vector3 _gpsAttachmentShift;

	// constants
	Eigen::Matrix<float, 3, 3> O33;
	Eigen::Matrix<float, 3, 4> O34;
	Eigen::Matrix<float, 4, 3> O43;
	Eigen::Matrix<float, 3, 3> E33;
};

