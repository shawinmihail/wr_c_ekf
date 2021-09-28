#include <iostream>
#include <stdio.h>
#include <vector>
#include <fstream>
#include <string.h>
#include <Eigen/Dense>
#include "EkfUnicorn.h"

using namespace Eigen;

int main(int argc,char** argv)
{
    EkfUnicornState Q0_diag;
    Q0_diag << 1e-6, 1e-6, 1e-6,  1e-4, 1e-4, 1e-4,  1e-7, 1e-7, 1e-7;
    EkfUnicornState P0_diag;
    P0_diag << 0.5, 0.5, 0.5, 0.04, 0.04, 0.04, 1e-4, 1e-4, 5e-2;
    Vector3f R_r_diag;
    R_r_diag << 1e-4, 1e-4, 1e-4;
    Vector3f R_v_diag;
    R_v_diag << 1e-4, 1e-4, 1e-4;
    Vector3f R_e_diag;
    R_e_diag << 1e-3, 1e-3, 1e-3;
    EkfUnicornInitParams initParams = {Q0_diag, P0_diag, R_r_diag, R_v_diag, R_e_diag};
    
    Vector3f drImuMaster(0.28, 0.0, 0.1);
    Vector3f drImuTarget(0, 0, 0);
    float K_a_smoothed = 0.025;

    
    EkfUnicorn ekfUnicorn;
    ekfUnicorn.setInitialParams(initParams);
    ekfUnicorn.setSensorsGeometry(drImuMaster, drImuTarget);
    ekfUnicorn.setAccSmoothParam(K_a_smoothed);
    
    Vector3f a(0.3, 0.4, 10);
    Vector3f w(0.3, 0.4, -0.1);
    Vector3f rZ(1, 2, -3);
    Vector3f vZ(0.12, 0.14, -0.15);
    float dt = 0.02;
    
    ekfUnicorn.updateImu(a, w);
    ekfUnicorn.defineInitialStateNoYaw(rZ, vZ);
    /*
    for (int i = 1; i <= 100; i++)
    {
        ekfUnicorn.predict(dt);
    }
    */
    //ekfUnicorn.predictNoYaw(dt);
    //ekfUnicorn.correct(rZ, vZ);
    ekfUnicorn.correctNoYaw(rZ, vZ);
    
    ekfUnicorn.printEstState();
    ekfUnicorn.printEstCov();

    
    return 0;
}
