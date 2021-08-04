#ifndef PLANAR_FLOATING_BASE
#define PLANAR_FLOATING_BASE

#include"PlanarRobot.h"
#include"FootholdPlan.h"

template<typename T>
class PlanarFloatingBase:public RobotBase<T>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
// private:
public:
    VecM<T,4> _foothold;
    VecM<T,2> _contact_state;
    VecM<T,xsize_FB> _xdot;
    MatMN<T,xsize_FB,xsize_FB> _Ac;
    MatMN<T,xsize_FB,usize_FB> _Bc; 
    FootholdPlanner<T> *_foothold_planner = nullptr;

protected:
    using RobotBase<T>::dynamics;
    using RobotBase<T>::dynamics_par;
    using RobotBase<T>::resetmap;
    using RobotBase<T>::resetmap_par;
    using RobotBase<T>::plan_foothold;

public:
    PlanarFloatingBase(FootholdPlanner<T>* fplanner, T dt = 0.001):
    RobotBase<T>(dt,xsize_FB,usize_FB,ysize_FB,ModelType::RBD)
    {
        _xdot.setZero();
        _Ac.setZero();
        _Bc.setZero();
        _foothold.setZero();
        _contact_state.setZero();
        _foothold_planner = fplanner;
    }

    void dynamics(VecM<T,xsize_FB> &x, VecM<T,usize_FB> &u, VecM<T, xsize_FB> &x_next, VecM<T,ysize_FB>&y, int mode) override;

    void resetmap(VecM<T,xsize_FB> &x, VecM<T,xsize_FB> &x_next, VecM<T,ysize_FB>&y, int mode) override;
    // void resetmap(DVec<T> &x, DVec<T> &x_next, DVec<T> &y, int mode = 1);

    void dynamics_par(VecM<T,xsize_FB>& x, VecM<T,usize_FB> &u, MatMN<T,xsize_FB,xsize_FB> &A, MatMN<T,xsize_FB,usize_FB> &B, 
                              MatMN<T,ysize_FB, xsize_FB> &C, MatMN<T, ysize_FB, usize_FB> &D, int mode ) override;
   
    void resetmap_par(VecM<T, xsize_FB> &x, MatMN<T,xsize_FB, xsize_FB> &Px, int mode) override;
    // void resetmap_par(DVec<T> &x, DMat<T> &Px, int mode = 1);

    void plan_foothold(DVec<T> &x,T stance_time, int mode) override;

};




#endif //PLANAR_FLOATING_BASE