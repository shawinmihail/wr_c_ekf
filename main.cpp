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

void main()
{

    SREKF srEkf;

    //std::string estDataPath("resources/est_states.floats");
    //MatlabDataParser estDater(estDataPath);
    //Vector16 estState;
    //float estTime = 0;

    //std::string actDataPath("resources/act_states.floats");
    //MatlabDataParser actDater(actDataPath);
    //Vector16 actState;
    //float actTime = 0;

    std::string mesDataPath("resources/mes_states.floats");
    MatlabDataParser mesDater(mesDataPath);
    Vector16 mesState;
    float mesTime = 0;

    int i = 0;
    while (true) {
        bool res = mesDater.next(mesState, mesTime);
        if (!res)
        {
            return;
        }

        if (i > 100)
        {
            return;
        }
        std::cout << i << "\n\n";

        srEkf.predictImu(mesState.segment(6, 3), mesState.segment(13, 3), 1e-3f);
        srEkf.correctPv(mesState.segment(0, 6));
        std::cout << srEkf.getEstState().transpose() << "\n\n";

        i++;
    }

    return;


    // logging init
    std::ofstream act_state_log;
    act_state_log.open("act_state_log.csv");

    std::ofstream mes_state_log;
    mes_state_log.open("mes_state_log.csv");

    std::ofstream est_state_log;
    est_state_log.open("est_state_log.csv");

    // init
    WheelRobotModel wheelRobotModel;
    MeasurmentModel mesModel;
    float dt = 1e-2f;
    Vector2 ctrl(0.5f, 0.1f);

    // loop
    for (int i = 0; i < 2000; i++)
    {
        // control input
        ctrl[0]=  0.90f; // v
        ctrl[1] = 0.25f; // u

        // moution
        wheelRobotModel.integrateMoution(ctrl, dt);
        //Vector6 state = wheelRobotModel.getState();
        Vector15 actState = wheelRobotModel.getFullState();

        // mesuarments
        mesModel.setData(actState);
        Vector15 mesState = mesModel.getMesState(); // r v a qv w
        
        // logging
        act_state_log << csvStrVect15(actState);
        mes_state_log << csvStrVect15(mesState);
        //est_state_log << csvStrVect9(estState);

        std::cout << i << std::endl;
    }

    act_state_log.close();
    mes_state_log.close();
    est_state_log.close();
}
