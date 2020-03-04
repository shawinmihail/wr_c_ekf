#include "MatlabDataParser.h"

#include <fstream>
#include <iostream>
#include <Eigen/Dense>


MatlabDataParser::MatlabDataParser(const std::string& path) :
    _length(0),
    _width(0),
    _data(NULL),
    _currentLineNum(0)
{
    std::ifstream ifs(path, std::ios::binary);

    ifs.read(reinterpret_cast<char*>(&_width), sizeof(int));
    ifs.read(reinterpret_cast<char*>(&_length), sizeof(int));

    const int datacount = (_length) * (_width);
    _data = (float*)malloc(sizeof(float) * datacount);
    if (!_data) {
        perror("Error allocating memory");
        abort();
    }

    int i = 0;
    while (ifs.read(reinterpret_cast<char*>(&_data[i]), sizeof(float))) {
        //std::cout << _data[i] << '\n';
        i++;
    }

}

bool MatlabDataParser::next(Vector16& state, float& dt)
{
    //Eigen::Map<Vector16> mp(_data + 17);
    //mp.

    if (_currentLineNum == _length - 1)
    {
        return false;
    }

    memcpy(&dt, _data + _currentLineNum * _width, sizeof(float));
    memcpy(&state(0), _data + 1 + _currentLineNum * _width, (_width - 1) * sizeof(float));
    _currentLineNum++;

    return true;
}

