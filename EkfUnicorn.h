#pragma once
#include <Eigen/Dense>
#include "Utils.h"

using namespace Eigen;
using namespace std;


// state X = [r v e]
#define EKF_UNICORN_DIM 9
typedef Matrix<float, EKF_UNICORN_DIM, 1> EkfUnicornState;
namespace EkfUnicornConstants
{
    Matrix<float, 3, 3> O33 (Matrix<float, 3, 3>::Zero());
	Matrix<float, 3, 4> O34 (Matrix<float, 3, 4>::Zero());
	Matrix<float, 3, 3> E33 (Matrix<float, 3, 3>::Identity());
};

using namespace EkfUnicornConstants;

struct EkfUnicornInitParams
{
    EkfUnicornState Q0_diag;
    EkfUnicornState P0_diag;
    Vector3f R_r_diag;
    Vector3f R_v_diag;
    Vector3f R_e_diag;
};


class EkfUnicorn
{
public:
	EkfUnicorn()
    {
        reinit();
    };
    
    void reinit()
    {
        // reset only state, for change params 
        // use setInitialParams(...) setSensorsGeometry(...) setAccSmoothParam(...)
        _a_smoothed = Vector3f::Zero();
        _a_smoothed_inited = false;
        _w = Vector3f::Zero();
        _X = EkfUnicornState::Zero();
    }
    
    void setInitialParams(const EkfUnicornInitParams& params)
    {
        _Q = params.Q0_diag.asDiagonal();
        _P = params.P0_diag.asDiagonal();
        _R_r = params.R_r_diag.asDiagonal();
        _R_v = params.R_v_diag.asDiagonal();
        _R_e = params.R_e_diag.asDiagonal();
        Vector9 R_diag;
        R_diag << params.R_r_diag, params.R_v_diag, params.R_e_diag;
        _R = R_diag.asDiagonal();
    };
    
    void setSensorsGeometry(const Vector3f& drImuMaster, const Vector3f& drImuTarget)
    {
        _drImuMaster = drImuMaster;
        _drImuTarget = drImuTarget;
        _drTargetMaster = _drImuMaster - _drImuTarget;
    };
    
    void setAccSmoothParam(float K_a_smoothed)
    {
        _K_a_smoothed = K_a_smoothed;
    };
    
    EkfUnicornState getEstState() const
	{
        return _X;
    };
    
	EkfUnicornState getEstTargetState() const
	{
        Vector3f r = _X.segment(0, 3);
        Vector3f v = _X.segment(3, 3);
        Vector3f e = _X.segment(6, 3);
        Vector4f q = quatFromEul(e);
        Vector3f rTarget = r + quatRotate(q, _drImuTarget);
        Vector3f vTarget = v + quatRotate(q, _drImuTarget.cross(_w));
        EkfUnicornState targetState;
        targetState << rTarget, vTarget, q;
        return targetState;
    };
    
    void printEstState() const
	{
        std::cout << "X: " << _X.transpose() << std::endl;
    };
    
    void printEstCov() const
	{
        std::cout << "D: " << _P.determinant() << std::endl;
        for (int i = 0; i < EKF_UNICORN_DIM; i++)
        {
            auto column = _P.col(i);
            std::cout << "column norm " << i << ": " << column.norm() << std::endl;
        }
            
        
        std::cout << "P:\n" << _P << std::endl;
    };
    
    void updateImu(const Vector3f& a, const Vector3f& w) // use corrected with calib data a, w
    {
        if (!_a_smoothed_inited)
        {
            _a_smoothed = a;
            _a_smoothed_inited = true;
        }
            
        _a_smoothed = smooth(a, _a_smoothed, _K_a_smoothed);
	    _w = w;
    }
    
    void defineInitialStateNoYaw(const Vector3f& r, const Vector3f& v)
    {
        // we use smoothed acc to define initial roll and pitch
        float roll = 0; float pitch = 0;
        rollPitchWithMeasuredAcc(_a_smoothed, roll, pitch);
        _X << r, v, roll, pitch, 0;
    };
    
	void predict(float dt)
    {
        /* state */
        Vector3f r0 = _X.segment(0, 3);
        Vector3f v0 = _X.segment(3, 3);
        Vector3f e0 = _X.segment(6, 3);
        Vector4f q0 = quatFromEul(e0);

        /* model */
        Vector3f r1 = r0 + v0 * dt;
        Vector3f v1 = v0 + quatRotate(q0, _w).cross(v0) * dt;
        Vector3f e1 = e0 + eulKinematics(e0, _w) * dt;
        _X << r1, v1, e1;

        /* covs */
        Matrix<float, 3, 3> Mvv = crossOperator(eulToMatrix(e0) * _w);
        Matrix<float, 3, 3> Mve = -crossOperator(v0) * eulRotLinearization(e0, _w);
        Matrix<float, 3, 3> Mee = eulKinLinearization(e0, _w);

        Matrix<float, 3, EKF_UNICORN_DIM> Fr;
        Matrix<float, 3, EKF_UNICORN_DIM> Fv;
        Matrix<float, 3, EKF_UNICORN_DIM> Fe;

        Fr << O33, E33, O33;
        Fv << O33, Mvv, Mve;
        Fe << O33, O33, Mee;

        Matrix<float, EKF_UNICORN_DIM, EKF_UNICORN_DIM> F;
        F << Fr, Fv, Fe;

        Matrix<float, EKF_UNICORN_DIM, EKF_UNICORN_DIM> I(
            Eigen::Matrix<float, EKF_UNICORN_DIM, EKF_UNICORN_DIM>::Identity());
        
        Matrix<float, EKF_UNICORN_DIM, EKF_UNICORN_DIM> PHI = I + F * dt;
        _P = PHI * _P * PHI.transpose() + _Q*dt;
    };
    
    void predictNoYaw(float dt)
    {
        /* 
        * Use predictNoYaw if yaw is not inited yet.
        * We can update pitch and roll with eulKinematics,
        * cause that derivateves do not depend on yaw.
        */
        
        /* state */
        Vector3f r0 = _X.segment(0, 3);
        Vector3f v0 = _X.segment(3, 3);
        Vector3f e0 = _X.segment(6, 3);
        Vector4f q0 = quatFromEul(e0);

        /* model */
        Vector3f r1 = r0 + v0 * dt;
        Vector3f v1 = v0 + quatRotate(q0, _w).cross(v0) * dt;
        Vector3f de = eulKinematics(e0, _w) * dt;
        de(2) = 0;
        Vector3f e1 = e0 + de;
        _X << r1, v1, e1;

        /* covs */
        Matrix<float, 3, 3> Mvv = crossOperator(eulToMatrix(e0) * _w);
        Matrix<float, 3, 3> Mve = -crossOperator(v0) * eulRotLinearization(e0, _w);
        Matrix<float, 3, 3> Mee = eulKinLinearization(e0, _w);

        Matrix<float, 3, EKF_UNICORN_DIM> Fr;
        Matrix<float, 3, EKF_UNICORN_DIM> Fv;
        Matrix<float, 3, EKF_UNICORN_DIM> Fe;

        Fr << O33, E33, O33;
        Fv << O33, Mvv, Mve;
        Fe << O33, O33, Mee;

        Matrix<float, EKF_UNICORN_DIM, EKF_UNICORN_DIM> F;
        F << Fr, Fv, Fe;
        
        Matrix<float, EKF_UNICORN_DIM-1, EKF_UNICORN_DIM-1> F_no_yaw = F.block<8, 8>(0,0);
        Matrix<float, EKF_UNICORN_DIM-1, EKF_UNICORN_DIM-1> P_no_yaw = _P.block<8, 8>(0,0);
        Matrix<float, EKF_UNICORN_DIM-1, EKF_UNICORN_DIM-1> Q_no_yaw = _Q.block<8, 8>(0,0);
        Matrix<float, EKF_UNICORN_DIM-1, EKF_UNICORN_DIM-1> I_no_yaw;
        I_no_yaw << Matrix<float, EKF_UNICORN_DIM-1, EKF_UNICORN_DIM-1>::Identity();
        Matrix<float, EKF_UNICORN_DIM-1, EKF_UNICORN_DIM-1> PHI_no_yaw = I_no_yaw + F_no_yaw * dt;
        P_no_yaw = PHI_no_yaw * P_no_yaw * PHI_no_yaw.transpose() + Q_no_yaw*dt;
        
        _P.block<8, 8>(0,0) = P_no_yaw;
    };
    
    void correct(const Vector3f& rZ, const Vector3f& vZ)
    {
        float roll_a = 0; float pitch_a = 0;
        rollPitchWithMeasuredAcc(_a_smoothed, roll_a, pitch_a);
        
        float yaw0 = _X(8);
        float pitch_v = 0; float yaw_v = 0;
        pitchYawWithMeasuredVel(vZ, yaw0, pitch_v, yaw_v);
        
        Vector3f eZ(roll_a, pitch_v, yaw_v);
        correct(rZ, vZ, eZ);
    }
    
    void correctNoYaw(const Vector3f& rZ, const Vector3f& vZ)
    {
        float roll_a = 0; float pitch_a = 0;
        rollPitchWithMeasuredAcc(_a_smoothed, roll_a, pitch_a);
        
        Vector3f eZ(roll_a, pitch_a, 0);
        std::cout << "eZ:\n" << eZ  << std::endl;
        correctNoYaw(rZ, vZ, eZ);
    }
    
    
private:
    
	void correct(const Vector3f& rZ, const Vector3f& vZ, const Vector3f& eZ)
    {
        Vector3f rX = _X.segment(0, 3);
        Vector3f vX = _X.segment(3, 3);
        Vector3f eX = _X.segment(6, 3);
        Vector4f qX = quatFromEul(eX);
        Vector9 Z;
        Z << rZ, vZ, eZ;

        // mes model
        // Z[rgnns vgnns e_mes]
        Vector3f Zr = rX + quatRotate(qX, _drImuMaster);
        Vector3f Zv = vX + quatRotate(qX, _w.cross(_drImuMaster));
        Vector3f Ze = eX;
        Vector9 Zx;
        Zx << Zr, Zv, Ze;
        Vector9 dz = Z - Zx;

        // H
        Matrix<float, 3, 3> Mre = eulRotLinearization(eX, _drImuMaster);
        Matrix<float, 3, 3> Mve = eulRotLinearization(eX, _w.cross(_drImuMaster));

        Matrix<float, 3, EKF_UNICORN_DIM> H1;
        Matrix<float, 3, EKF_UNICORN_DIM> H2;
        Matrix<float, 3, EKF_UNICORN_DIM> H3;
        Matrix<float, 9, EKF_UNICORN_DIM> H;
        H1 << E33, O33, Mre;
        H2 << O33, E33, Mve;
        H3 << O33, O33, E33;
        H << H1, H2, H3;

        // ordinary KF
        Matrix<float, 9, 9> Rk = _R + H * _P * H.transpose();
        Matrix<float, EKF_UNICORN_DIM, 9> K = _P * H.transpose() * (Rk.inverse());
        Matrix<float, EKF_UNICORN_DIM, EKF_UNICORN_DIM> dP = K * H * _P;
        _P = _P - dP;
        EkfUnicornState dx = K * dz;
        _X = _X + dx;
    };
    
	void correctNoYaw(const Vector3f& rZ, const Vector3f& vZ, Vector3f eZ)
    {
        /* 
        * Use correctNoYaw if yaw is not inited yet
        * to update pos vel with gnns 
        * and roll pitch with imu
        * 
        * suppose eZ(2) = 0;
        */
        
        Vector3f rX = _X.segment(0, 3);
        Vector3f vX = _X.segment(3, 3);
        Vector3f eX = _X.segment(6, 3);
        Vector4f qX = quatFromEul(eX);
        Vector9 Z;
        eZ(2) = 0;
        Z << rZ, vZ, eZ;

        // mes model
        // Z[rgnns vgnns e_mes]
        Vector3f Zr = rX + quatRotate(qX, _drImuMaster);
        Vector3f Zv = vX + quatRotate(qX, _w.cross(_drImuMaster));
        Vector3f Ze = eX;
        Vector9 Zx;
        Zx << Zr, Zv, Ze;
        Vector9 dz = Z - Zx;
        Vector8 dz_no_yaw = dz.segment(0, 8);

        // H
        Matrix<float, 3, 3> Mre = eulRotLinearization(eX, _drImuMaster);
        Matrix<float, 3, 3> Mve = eulRotLinearization(eX, _w.cross(_drImuMaster));

        Matrix<float, 3, EKF_UNICORN_DIM> H1;
        Matrix<float, 3, EKF_UNICORN_DIM> H2;
        Matrix<float, 3, EKF_UNICORN_DIM> H3;
        Matrix<float, 9, EKF_UNICORN_DIM> H;
        H1 << E33, O33, Mre;
        H2 << O33, E33, Mve;
        H3 << O33, O33, E33;
        H << H1, H2, H3;
        Matrix<float, 8, EKF_UNICORN_DIM-1> H_no_yaw = H.block<8, EKF_UNICORN_DIM-1>(0, 0);

        // ordinary KF
        Matrix<float, 8, 8> R_no_yaw = _R.block<8, 8>(0, 0);
        Matrix<float, 8, 8> P_no_yaw = _P.block<8, 8>(0, 0);
        
        Matrix<float, 8, 8> Rk_no_yaw = R_no_yaw + H_no_yaw * P_no_yaw * H_no_yaw.transpose();
        Matrix<float, EKF_UNICORN_DIM-1, 8> K_no_yaw = P_no_yaw * H_no_yaw.transpose() * (Rk_no_yaw.inverse());
        Matrix<float, EKF_UNICORN_DIM-1, EKF_UNICORN_DIM-1> dP_no_yaw = K_no_yaw * H_no_yaw * P_no_yaw;
        P_no_yaw = P_no_yaw - dP_no_yaw;
        EkfUnicornState dx;
        dx << K_no_yaw * dz_no_yaw, 0;
        _X = _X + dx;
        _P.block<8, 8>(0,0) = P_no_yaw;
    };

private:
	Vector3f smooth(const Vector3& sample, const Vector3& smoothed, float K)
    {
        Vector3f res = smoothed + K * (sample - smoothed);
        return res;
    };
    

    
    void rollPitchWithMeasuredAcc(const Vector3f& acc, float& roll, float& pitch)
    {
        // we get bow with measured acc g dir
        
        roll = 0;
        pitch = 0;
        
        float n_acc = acc.norm();
        if (n_acc < 1)
        {
            return;
        }
        
        Vector3f dir = acc / n_acc;
        pitch = asin(-dir[0]);
        roll = asin(dir[1] / cos(pitch));
        return;
    };
    
    void pitchYawWithMeasuredVel(const Vector3f& v, float yaw0, float& pitch, float& yaw)
    {
        // we use velocity to define yaw and pitch
        // we use incremental yaw (-inf inf) to avoid correction gap in +- pi
                
        Vector3f ex(1,0,0);
        Vector4f q = quatBetweenVectors(ex, v);
        Vector3f rpy = quatToEul(q);
        pitch = rpy[1];
        float yawCutPmPi = rpy[2];
        
        float dyaw = shortestRotation(yaw0 ,yawCutPmPi);
        yaw = yaw0 + dyaw;
        
        return;
    };


private:
	EkfUnicornState _X;
    
	Matrix<float, EKF_UNICORN_DIM, EKF_UNICORN_DIM> _P;
	Matrix<float, EKF_UNICORN_DIM, EKF_UNICORN_DIM> _Q;
	Matrix<float, 3, 3> _R_r;
	Matrix<float, 3, 3> _R_v;
	Matrix<float, 3, 3> _R_e;
	Matrix<float, 9, 9> _R;

	Vector3f _drImuMaster;
	Vector3f _drImuTarget;
	Vector3f _drTargetMaster;
    
    Vector3f _a_smoothed;
    bool _a_smoothed_inited;
	float _K_a_smoothed;
	Vector3f _w;

EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

