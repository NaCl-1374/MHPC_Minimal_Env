#ifndef PLANAR_ROBOT_H
#define PLANAR_ROBOT_H

#include "MHPC_CPPTypes.h"

enum class ModelType
{
    FULL,       // Whole body 
    RBD,        // rigid body
    CENTROIDAL  // potato dyanmics
};
const size_t xsize_WB = 14, usize_WB = 4, ysize_WB = 4, qsize_WB=7;
const size_t xsize_FB = 6, usize_FB = 4, ysize_FB = 4, qsize_FB=3;

template<typename T>
class RobotBase
{
public:
    /* data */
    T _dt;
    const size_t _qsize, _xsize, _usize, _ysize;   // dof, state, control, output dimensions
    ModelType _mtype;

public:
    RobotBase(T dt, size_t xsize, size_t usize, size_t ysize, ModelType mtype): _dt(dt),
    _qsize(xsize/2),_xsize(xsize), _usize(usize), _ysize(ysize) {_mtype = mtype;}

    virtual ~RobotBase() = default;

    ModelType get_model_type() {return _mtype;}
    
    virtual void dynamics(VecM<T,xsize_WB> &x, VecM<T,usize_WB> &u, VecM<T, xsize_WB> &x_next, VecM<T,ysize_WB>&y, int mode = 1) {}

    virtual void dynamics(VecM<T,xsize_FB> &x, VecM<T,usize_FB> &u, VecM<T, xsize_FB> &x_next, VecM<T,ysize_FB>&y, int mode = 1) {}

    virtual void resetmap(VecM<T,xsize_WB> &x, VecM<T,xsize_WB> &x_next, VecM<T,ysize_WB>&y, int mode = 1) {}

    virtual void resetmap(VecM<T,xsize_FB> &x, VecM<T,xsize_FB> &x_next, VecM<T,ysize_FB>&y, int mode = 1) {}

    // virtual void resetmap(DVec<T> &x, DVec<T> &x_next, DVec<T>&y, int mode = 1) {}

    virtual void dynamics_par(VecM<T,xsize_WB>& x, VecM<T,usize_WB> &u, MatMN<T,xsize_WB,xsize_WB> &A, MatMN<T,xsize_WB,usize_WB> &B, 
                              MatMN<T,ysize_WB, xsize_WB> &C, MatMN<T, ysize_WB, usize_WB> &D, int mode = 1) {}

    virtual void dynamics_par(VecM<T,xsize_FB>& x, VecM<T,usize_FB> &u, MatMN<T,xsize_FB,xsize_FB> &A, MatMN<T,xsize_FB,usize_FB> &B, 
                              MatMN<T,ysize_FB, xsize_FB> &C, MatMN<T, ysize_FB, usize_FB> &D, int mode = 1) {}
   
    virtual void resetmap_par(VecM<T, xsize_WB> &x, MatMN<T,xsize_WB, xsize_WB> &Px, int mode = 1) {}

    virtual void resetmap_par(VecM<T, xsize_FB> &x, MatMN<T,xsize_FB, xsize_FB> &Px, int mode = 1) {}
    // virtual void resetmap_par(DVec<T> &x, DMat<T> &Px, int mode = 1) {}

    virtual void plan_foothold(VecM<T, xsize_WB> &x,T stance_time, int mode = 1){}

    virtual void plan_foothold(DVec<T> &x,T stance_time, int mode = 1){}

};


#endif //PLANAR_ROBOT_H