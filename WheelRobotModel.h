#pragma once
#include "Definitions.h"

class WheelRobotModel
{
public:
	WheelRobotModel();
	void integrateMoution(const Vector2& ctrl, float dt);
	Vector6 getState();
	Vector15 getFullState();

private:
	void processCtrlInput(const Vector2& ctrl, float dt);
	void updateStateDot(float dt);
	void updateState(float dt);
	void updateFullState();
private:
	Vector15 _fullState; // r3 v3 a3 qv3 w3
	Vector6 _state; // x y th x_dot y_dot th_dot
	Vector6 _stateDot;
	Vector6 _stateDotLast;
	Vector2 _ctrl; // v u
};

