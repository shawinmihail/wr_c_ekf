#pragma once
#include <string>
#include "Definitions.h"

class MatlabDataParser {
public:
	MatlabDataParser(const std::string& path);
	bool next(Vector22& state, float& dt);
private:
	int _length;
	int _width;
	float* _data;

	int _currentLineNum;
};