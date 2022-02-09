#pragma once
#include "EkfWgsMath.h"

class EKF4_WGS
{
public:
	EKF4_WGS();
    void initParams();
    bool initState(const Eigen::Vector3d& r0, const Eigen::Vector3d& s1, const Eigen::Vector3d& s2, double slaves_accuracy);
	void predict(double dt);
	void correctRV(const Eigen::Vector3d& r, const Eigen::Vector3d& v);
	void correctQ2(const Eigen::Vector3d& dr1, const Eigen::Vector3d& dr2);
	void setImu(const Eigen::Vector3d& a, const Eigen::Vector3d& w);
    void setQImuCalib(const Eigen::Vector4d& q);
    void setBiasesImuCalib(const Eigen::Vector3d& da, const Eigen::Vector3d& dw);
    void setSlavesCalib(const Eigen::Vector3d& slave1, const Eigen::Vector3d& slave2);
    void setSensorsGeometry(const Eigen::Vector3d& drImuMaster, const Eigen::Vector3d& drImuTarget);
    //
	Eigen::Matrix<double, 10, 1> getEstState();
	Eigen::Matrix<double, 10, 1> getEstTargetState();
    Eigen::Vector3d getDrSlave1();
    Eigen::Vector3d getDrSlave2();
    void getImuCorrected(Eigen::Vector3d& a, Eigen::Vector3d& w);
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

private:
	Eigen::Vector3d smooth(const Eigen::Vector3d& sample, const Eigen::Vector3d& smoothed, double K);
	Eigen::Vector3d _a_smoothed;
	Eigen::Vector3d _w_smoothed;
	double _K_a_smoothed;
	double _K_w_smoothed;
    Eigen::Vector4d _qImuCalib;
    Eigen::Vector3d _daImuCalib;
    Eigen::Vector3d _dwImuCalib;
    
public:
    bool attitudeWithBaseLineNumeric(const Eigen::Vector3d& s1, const Eigen::Vector3d& s2, double slaves_accuracy, Eigen::Vector4d& q_out);
    
private:
	Eigen::Matrix<double, 10, 1> _X;
	Eigen::Matrix<double, 10, 10> _Q;
	Eigen::Matrix<double, 10, 10> _P;
	Eigen::Matrix<double, 6, 6> _R_RV;
	Eigen::Matrix<double, 6, 6> _R_Q2;

	Eigen::Vector3d _drImuMaster;
	Eigen::Vector3d _drImuTarget;
	Eigen::Vector3d _drTargetMaster;
	Eigen::Vector3d _drSlave1;
	Eigen::Vector3d _drSlave2;

	// constants
	Eigen::Matrix<double, 3, 3> O33;
	Eigen::Matrix<double, 3, 4> O34;
	Eigen::Matrix<double, 4, 3> O43;
	Eigen::Matrix<double, 3, 3> E33;
};

