#pragma once

#include "Utils.h"
#include <random>
#include "Definitions.h"

class MeasurmentModel
{
public:
	MeasurmentModel();
	void setData(const Vector15& fullState);
	Vector15 getMesState();
private:
	Vector15 _mesState; // r v a qv w

	std::default_random_engine _generator;
	std::normal_distribution<float> _imuAccDist;
	std::normal_distribution<float> _imuRotVelDist;
	std::normal_distribution<float> _gpsPosDist;
	std::normal_distribution<float> _gpsVelDist;
};

