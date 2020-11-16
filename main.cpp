#pragma once

#include <iostream>
#include <fstream>
#include <sstream>

#include "Definitions.h"
#include "Utils.h"
#include "Ekf4.h"

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
    EKF4 ekf;

    // read mes
    std::string mesDataPath("resources/mes_states.floats");
    MatlabDataParser mesDater(mesDataPath);
    MatlabDataParserVector mesState;
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

        if (i > 999)
        {
            return;
        }

        std::cout << i << "\n\n";

        ekf.setImu(mesState.segment(12, 3), mesState.segment(15, 3));
        ekf.predict(1e-1f);
        ekf.correctRV(mesState.segment(0, 6));
        //ekf.correctU(mesState.segment(3, 3));
        //srEkf.correctA(mesState.segment(6, 3));
        ekf.correctQ2(mesState.segment(6, 3), mesState.segment(9, 3));
        //std::cout << mesState.segment(16, 3) << "\n\n";
        //std::cout << ekf.getEstState().segment(6, 4).transpose() << "\n\n";
        //std::cout << mesState.transpose() << "\n\n";

        /* logs */
        Ekf4_fullState X = ekf.getEstState();
        Vector3 r = X.segment(0, 3);
        Vector3 v = X.segment(3, 3);
        Vector4 q = X.segment(6, 4);
        Vector3 drSlave1;
        drSlave1 << 0.73f, 0.23f, 0.0f;
        Vector3 drSlave2;
        drSlave2 << 0.73f, -0.23f, 0.0f;
        Vector3 drSlave1Est = quatRotate(q, drSlave1);
        Vector3 drSlave2Est = quatRotate(q, drSlave2);
        Eigen::Matrix<float, 24, 1> logMsg;
        logMsg << r, v, drSlave1Est, drSlave2Est, mesState.segment(0, 6), mesState.segment(6, 6);
        est_state_log << csvStrVect<24>(logMsg);

        i++;
    }

    // postproc
    est_state_log.close();
}
