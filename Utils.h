#pragma once
#include <stdio.h>
#include "Definitions.h"

// quat

Vector4 quatFromEul(const Vector3& eul);

Vector3 quatToEul(const Vector4& q);

Vector4 quatMultiply(const Vector4& q1, const Vector4& q2);

Vector4 quatInverse(const Vector4& q);

Vector3 quatRotate(const Vector4& q, const Vector3& v);

Eigen::Matrix<float, 3, 3> quatToMatrix(const Vector4& q);

Eigen::Matrix<float, 3, 3> crossOperator( const Vector3& v);

Vector4 quatBetweenVectors(const Vector3& v1, const Vector3& v2);

float shortestRotation(float from, float to);

// lin

Eigen::Matrix<float, 3, 4> quatRotateLinearizationQ(const Vector4& q, const Vector3& v);

Eigen::Matrix<float, 4, 4> poissonEqLinearizationQ(const Vector3& w);

Eigen::Matrix<float, 4, 3> poissonEqLinearizationW(const Vector4& q, const Vector3& w);

Eigen::Matrix<float, 1, 3> normVect3Linearization(const Vector3& v);


Eigen::Matrix<float, 3, 3> eulToMatrix(const Vector3& e);

Vector3 eulKinematics(const Vector3& e, const Vector3& w);

Eigen::Matrix<float, 3, 3> eulRotLinearization(const Vector3& e, const Vector3& w);

Eigen::Matrix<float, 3, 3> eulKinLinearization(const Vector3& e, const Vector3& w);



