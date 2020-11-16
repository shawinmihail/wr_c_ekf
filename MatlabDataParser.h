#pragma once
#include <string>
#include "Definitions.h"

typedef Eigen::Matrix<float, 18, 1> MatlabDataParserVector;
class MatlabDataParser
{
public:
	MatlabDataParser(const std::string& path);
	bool next(MatlabDataParserVector& state, float& dt);
private:
	int _length;
	int _width;
	float* _data;

	int _currentLineNum;
};