#include "WheelRobotModel.h"
#include "Utils.h"
#include <math.h>

WheelRobotModel::WheelRobotModel() :
_state(Vector6::Zero()),
_stateDot(Vector6::Zero()),
_stateDotLast(Vector6::Zero()),
_fullState(Vector15::Zero()),
_ctrl(Vector2::Zero())
{

}

void WheelRobotModel::updateFullState()
{

	// r
	_fullState[0] = _state[0]; // x
	_fullState[1] = _state[1]; // y
	_fullState[2] = 0.0f;      // z
	// v
	_fullState[3] = _state[3]; // vx
	_fullState[4] = _state[4]; // vy
	_fullState[5] = 0.0f;      // vz
	// a
	_fullState[6] = _stateDot[3]; // ax
	_fullState[7] = _stateDot[4]; // ay
	_fullState[8] = 0.0f;         // az
	// qv
	Vector3 eul(0.f, 0.f, _state(2, 0));  // eul = [0;0;th];
	Vector4 q = quatFromEul(eul);
	Vector3 qv = quatToQuatVec(q);
	_fullState(9, 0) =  qv[0];  // qx
	_fullState(10, 0) = qv[1];  // qy
	_fullState(11, 0) = qv[2];  // qz
	// w = [0;0;th_dot]
	_fullState[12] = 0.0f;         // wx
	_fullState[13] = 0.0f;         // wy
	_fullState[14] = _stateDot[2]; // wz
}

void WheelRobotModel::integrateMoution(const Vector2& ctrl, float dt)
{
	processCtrlInput(ctrl, dt);
	updateState(dt);
	updateStateDot(dt);
	//printVect6(_stateDot);
	updateFullState();
}

void WheelRobotModel::processCtrlInput(const Vector2& ctrl, float dt)
{
	float Kv = 1.0f;
	float Ku = 1.0f;
	float dv = ctrl[0] - _ctrl[0];
	float du = ctrl[1] - _ctrl[1];

	_ctrl[0] = _ctrl[0] + Kv * dv * dt;
	_ctrl[1] = _ctrl[1] + Ku * du * dt;
}

void WheelRobotModel::updateState(float dt)
{
	_state = _state + _stateDot * dt;
}

void WheelRobotModel::updateStateDot(float dt)
{
	_stateDot[0] = _ctrl[0] * cos(_state[2]); // x_dot = v cos (th)
	_stateDot[1] = _ctrl[0] * sin(_state[2]); // y_dot = v sin (th)
	_stateDot[2] = _ctrl[0] * _ctrl[1];       // th_dot = v u

	// numerical
	_stateDot[3] = (_stateDot[0] - _stateDotLast[0]) / dt; // (xdot_k - xdot_k-1) / dt
	_stateDot[4] = (_stateDot[1] - _stateDotLast[1]) / dt;
	_stateDot[5] = (_stateDot[2] - _stateDotLast[2]) / dt;

	_stateDotLast = _stateDot;
}

Vector6 WheelRobotModel::getState()
{
	return _state;
}

Vector15 WheelRobotModel::getFullState()
{
	return _fullState;
}