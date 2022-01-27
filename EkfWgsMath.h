#pragma once
#include <stdio.h>
#include <Eigen/Dense>

// quat

Eigen::Vector4d quatFromEul(const Eigen::Vector3d& eul);

Eigen::Vector3d quatToEul(const Eigen::Vector4d& q);

Eigen::Vector4d quatMultiply(const Eigen::Vector4d& q1, const Eigen::Vector4d& q2);

Eigen::Vector4d quatInverse(const Eigen::Vector4d& q);

Eigen::Vector3d quatRotate(const Eigen::Vector4d& q, const Eigen::Vector3d& v);

Eigen::Matrix<double, 3, 3> quatToMatrix(const Eigen::Vector4d& q);

Eigen::Matrix<double, 3, 3> crossOperator( const Eigen::Vector3d& v);

Eigen::Vector4d quatBetweenVectors(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2);

double shortestRotation(double from, double to);

// lin

Eigen::Matrix<double, 3, 4> quatRotateLinearizationQ(const Eigen::Vector4d& q, const Eigen::Vector3d& v);

Eigen::Matrix<double, 4, 4> poissonEqLinearizationQ(const Eigen::Vector3d& w);

Eigen::Matrix<double, 4, 3> poissonEqLinearizationW(const Eigen::Vector4d& q, const Eigen::Vector3d& w);

Eigen::Matrix<double, 1, 3> normVect3Linearization(const Eigen::Vector3d& v);

Eigen::Matrix<double, 3, 3> eulRotLinearization(const Eigen::Vector3d& e, const Eigen::Vector3d& w);

Eigen::Matrix<double, 3, 3> eulKinLinearization(const Eigen::Vector3d& e, const Eigen::Vector3d& w);

// eul

Eigen::Vector3d eulKinematics(const Eigen::Vector3d& e, const Eigen::Vector3d& w);

Eigen::Matrix<double, 3, 3> eulToMatrix(const Eigen::Vector3d& e);

// WGS

Eigen::Matrix<double, 3, 3> wgsToEnuRmat(double lat, double lon);

Eigen::Vector4d wgsToEnuQuat(double lat, double lon);





