#pragma once
#include "Definitions.h"
#include "Utils.h"

// state X = [r v q]
#define EKF4_STATE_DIM 9
typedef Eigen::Matrix<float, EKF4_STATE_DIM, 1> Ekf4_state;
typedef Eigen::Matrix<float, EKF4_STATE_DIM + 1, 1> Ekf4_fullState;


class EKF4
{
public:
	EKF4();
    void initParams();
    void reset(const Vector3& r0);
    void calibSlavesWithSample(const Vector3& dr1, const Vector3& dr2);
	void predict(float dt);
	void correctRV(const Vector6& rv);
	void correctRV(const Vector3& r, const Vector3& v);
	void correctU(const Vector3& v);
	void correctA(const Vector3& a);
	void correctQ2(const Vector3& dr1, const Vector3& dr2);
	void setImu(const Vector3& a, const Vector3& w);
	Ekf4_fullState getEstState();

private:
	Vector3 smooth(const Vector3& sample, const Vector3& smoothed, float K);
	Vector3 _a_smoothed;
	Vector3 _w_smoothed;
	float _K_a_smoothed;
	float _K_w_smoothed;

private:
	Ekf4_fullState _X;
	Eigen::Matrix<float, EKF4_STATE_DIM, EKF4_STATE_DIM> _Q;
	Eigen::Matrix<float, EKF4_STATE_DIM, EKF4_STATE_DIM> _P;

	Eigen::Matrix<float, 6, 6> _R_RV;
	Eigen::Matrix<float, 6, 6> _R_Q2;
	Eigen::Matrix<float, 3, 3> _R_U;
	Eigen::Matrix<float, 3, 3> _R_A;

	Vector3 _drImuGnns;
	Vector3 _drBshcGnns;
	Vector3 _drSlave1;
	Vector3 _drSlave2;

	// constants
	Eigen::Matrix<float, 3, 3> O33;
	Eigen::Matrix<float, 3, 4> O34;
	Eigen::Matrix<float, 4, 3> O43;
	Eigen::Matrix<float, 3, 3> E33;
};

