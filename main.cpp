#pragma once

#include <iostream>
#include <fstream>
#include <sstream>

#include "Definitions.h"
#include "Utils.h"
#include "WheelRobotModel.h"
#include "MeasurmentModel.h"
#include "SREKF.h"

#include "MatlabDataParser.h"

#include <Eigen/Dense>

// csv
template <int n>
std::string csvStrVect(const Eigen::Matrix<float, n, 1> &vect)
{
    std::ostringstream os;
    for (int i = 0; i < n; i++)
    {
        os << vect[i];
        os << ",";
    }
    os << "\n";
    std::string str(os.str());
    return str;
}

void main()
{
    // init
    std::ofstream est_state_log;
    est_state_log.open("src/est_state_log.csv");
    SREKF srEkf;

    // read mes
    std::string mesDataPath("resources/mes_states.floats");
    MatlabDataParser mesDater(mesDataPath);
    Vector16 mesState;
    float mesTime = 0;

    // test
    int i = 0;
    while (true)
    {
        bool res = mesDater.next(mesState, mesTime);
        if (!res)
        {
            return;
        }

        if (i > 19999)
        {
            return;
        }

        std::cout << i << "\n\n";

        srEkf.predictImu(mesState.segment(6, 3), mesState.segment(13, 3), 1e-3f);
        //srEkf.correctPv(mesState.segment(0, 6));
        srEkf.correctZ(mesState.segment(0, 3), mesState.segment(3, 3), mesState.segment(6, 3));
        //std::cout << srEkf.getEstState().transpose() << "\n\n";

        est_state_log << csvStrVect<16>(srEkf.getEstState());

        i++;
    }

    // postproc
    est_state_log.close();
}
