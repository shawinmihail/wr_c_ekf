#include "MeasurmentModel.h"

MeasurmentModel::MeasurmentModel() :
_mesState(Vector15::Zero()),

_generator(200),
_imuAccDist(0.0, 0.05f),
_imuRotVelDist(0.0f, 0.005f),
_gpsPosDist(0.0f, 0.5f),
_gpsVelDist(0.0f, 0.05f)
{

}

void MeasurmentModel::setData(const Vector15& fullState /* r v a qv w*/)
{
	// pose vel in local inertial frame

	_mesState[0] = fullState[0] + _gpsPosDist(_generator);
	_mesState[1] = fullState[1] + _gpsPosDist(_generator);
	_mesState[2] = fullState[2] + _gpsPosDist(_generator);

	_mesState[3] = fullState[3] + _gpsVelDist(_generator);
	_mesState[4] = fullState[4] + _gpsVelDist(_generator);
	_mesState[5] = fullState[5] + _gpsVelDist(_generator);

	// attitude
	Vector3 qv(fullState[9], fullState[10], fullState[11]);
	Vector4 q = quatVecToQuat(qv);
	_mesState[9] = fullState[9];
	_mesState[10] = fullState[10];
	_mesState[11] = fullState[11];

	// imu in body frame
	float gz = -10.0f;
	Vector3 a(fullState[6], fullState[7], fullState[8] - gz);
	Vector4 quatDual = quatInverse(q);
	Vector3 aB = quatRotate(quatDual, a);
	_mesState[6] = aB[0] + _imuAccDist(_generator);
	_mesState[7] = aB[1] + _imuAccDist(_generator);
	_mesState[8] = aB[2] + _imuAccDist(_generator);

	// suppose rot rate the same in I and B frames (only Z rotation)
	Vector3 w(3, 1);
	_mesState[12] = fullState[12] + _imuRotVelDist(_generator);
	_mesState[13] = fullState[13] + _imuRotVelDist(_generator);
	_mesState[14] = fullState[14] + _imuRotVelDist(_generator);
}

Vector15 MeasurmentModel::getMesState()
{
	return _mesState;
}