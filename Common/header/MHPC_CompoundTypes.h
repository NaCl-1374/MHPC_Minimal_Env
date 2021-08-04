#ifndef MHPC_COMPOUNDTYPES_H
#define MHPC_COMPOUNDTYPES_H

#include "MHPC_CPPTypes.h"
#include <chrono>

template<typename T, size_t xsize, size_t usize, size_t ysize>
struct ModelState
{   
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    VecM<T, xsize> x;
    VecM<T, usize> u;
    VecM<T, ysize> y;

    ModelState() {Zeros();}
    void Zeros()
    {
        x.setZero();
        u.setZero();
        y.setZero();
    }
};

template<typename T, size_t xsize, size_t usize, size_t ysize>
struct DynDerivative
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    MatMN<T, xsize, xsize> A;
    MatMN<T, xsize, usize> B;
    MatMN<T, ysize, xsize> C;
    MatMN<T, ysize, usize> D;

    DynDerivative() {Zeros();}
    void Zeros()
    {
        A.setZero();
        B.setZero();
        C.setZero();
        D.setZero();        
    }

};

template<typename T, size_t xsize, size_t usize, size_t ysize>
struct RCostStruct
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    T l;
    VecM<T, xsize> lx;
    VecM<T, usize> lu;
    VecM<T, ysize> ly;
    MatMN<T, xsize, xsize> lxx;
    MatMN<T, usize, xsize> lux;
    MatMN<T, usize, usize> luu;
    MatMN<T, ysize, ysize> lyy;

    RCostStruct(){Zeros();}
    void Zeros()
    {
        l = 0;
        lx.setZero();
        lu.setZero();
        ly.setZero();
        lxx.setZero();
        luu.setZero();
        lux.setZero();
        lyy.setZero();
    }
};

template<typename T, size_t xsize>
struct TCostStruct
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    T Phi;
    VecM<T, xsize> Phix;
    MatMN<T, xsize, xsize> Phixx;

    TCostStruct(){Zeros();}
    void Zeros()
    {
        Phi = 0;
        Phix.setZero();
        Phixx.setZero();
    }
};

template<typename T, size_t xsize, size_t usize>
struct CostToGoStruct
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    VecM<T, xsize> G; // Vx
    MatMN<T, xsize, xsize> H; // Vxx
    VecM<T, usize> du;
    MatMN<T, usize, xsize> K;
    VecM<T, xsize> Qx;
    VecM<T, usize> Qu;
    MatMN<T, xsize, xsize> Qxx;
    MatMN<T, usize, xsize> Qux;
    MatMN<T, usize, usize> Quu;

    CostToGoStruct(){Zeros();}
    void Zeros()
    {
        G.setZero();
        H.setZero();
        K.setZero();
        du.setZero();
        Qx.setZero();
        Qu.setZero();
        Qxx.setZero();
        Qux.setZero();
        Quu.setZero();
    }

    template<size_t ysize>
    void compute_Qfunction(RCostStruct<T, xsize, usize, ysize> &rcost, DynDerivative<T, xsize, usize, ysize> &dynpar, 
                            VecM<T, xsize> &Gnext, MatMN<T, xsize, xsize> &Hnext)
    {
        // compute Q function approximation
        Qx = rcost.lx + dynpar.A.transpose() * Gnext + dynpar.C.transpose() * rcost.ly;
        Qu = rcost.lu + dynpar.B.transpose() * Gnext + dynpar.D.transpose() * rcost.ly;
        Qxx = rcost.lxx + dynpar.C.transpose() * rcost.lyy * dynpar.C + dynpar.A.transpose() * Hnext * dynpar.A;
        Quu = rcost.luu + dynpar.D.transpose() * rcost.lyy * dynpar.D + dynpar.B.transpose() * Hnext * dynpar.B;
        Qux = rcost.lux + dynpar.D.transpose() * rcost.lyy * dynpar.C + dynpar.B.transpose() * Hnext * dynpar.A;
    }

    T valuefunction_update()
    {
        T dV; // expected cost change at current time step (not summed up)

        // Symmetrize Quu_inv and Qxx. Numerical issue would occur otherwise.
        MatMN<T,usize, usize> Quu_inv = (Quu.inverse() + Quu.inverse().transpose())/2;
        Qxx = (Qxx + Qxx.transpose())/2;

        // compute value function approximation
        du = -Quu_inv * Qu;
        K = -Quu_inv * Qux;
        G = Qx - Qux.transpose() * Quu_inv * Qu;
        H = Qxx - Qux.transpose() * Quu_inv * Qux;
        dV = -Qu.transpose() * Quu.inverse() * Qu;

        return dV;
    }
};

template<typename T, size_t xsize>
struct TConstrStruct
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    T h;  // put all terminal constraint in a vector
    VecM<T, xsize> hx;
    MatMN<T,xsize, xsize> hxx;
   
    T get_violation_normsquare()
    {
        return h*h;
    }

    TConstrStruct() {Zeros();}

    void Zeros()
    {
        h = 0;
        hx.setZero();
        hxx.setZero();
    }
};

template<typename T, size_t xsize, size_t usize, size_t ysize>
struct IneqConstrStruct
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // probably not needed in this structure
    T g;
    VecM<T, xsize> gx;
    VecM<T, usize> gu;
    VecM<T, ysize> gy;
    MatMN<T, xsize, xsize> gxx;
    MatMN<T, usize, usize> guu;
    MatMN<T, ysize, ysize> gyy;

    IneqConstrStruct() {Zeros();}
    void Zeros()
    {
        g = 0;
        gx.setZero();
        gu.setZero();
        gy.setZero();
        gxx.setZero();
        guu.setZero();
        gyy.setZero();
    }

};

template<typename T>
struct HSDDP_OPTION
{
    T alpha = 0.1;               // line search udpate paramer
    T gamma = 0.01;              // scale the expected cost reduction
    T update_penalty = 8;        // penalty update parameter
    T update_relax = 0.1;        // relaxation parameter udpate parameter
    T update_regularization = 2; // regularization parameter update parameter
    T update_ReB = 7;            // update barrier function weighting
    T max_DDP_iter = 3;          // maximum inner loop iteration
    T max_AL_iter = 2;           // maximum outer loop iteration
    T DDP_thresh = 1e-03;        // inner loop convergence threshhold
    T AL_thresh = 1e-03;         // outer loop convergence threshhold
    bool AL_active = 1;               // activate terminal constraint
    bool ReB_active = 1;              // activate path constraint
    bool smooth_active = 0;           // activate control smoothness penalization
};

template<typename T>
struct AL_REB_PARAMETER
{
    bool al_empty = true;      // assume no switching constraint
    bool reb_empty = true;     // assume no inequality constraint
    T sigma = 0;      // penalty parameter
    DVec<T> lambda;           // lagrangian parameters
    DVec<T> delta;            // relaxation parameters
    DVec<T> delta_min;        // lower bound on relaxation parameters
    DVec<T> eps_ReB;          // weighting parameters for the barrier function
    T eps_smooth = 1; // weighting parameter for control smoothness

    AL_REB_PARAMETER(){}; // default constructor

    AL_REB_PARAMETER(size_t num_eq, size_t num_ineq)
    {
        lambda.setZero(num_eq);
        delta = 0.1 * DVec<T>::Ones(num_ineq);
        delta_min = 0.01 * DVec<T>::Ones(num_ineq);
        eps_ReB.setZero(num_ineq);
    }
};

struct USRCMD
{
    float vel, height, roll, pitch, yaw;
};

struct  MHPCUserParameters
{
    int n_wbphase = 4;
    int n_fbphase = 4;
    float dt_wb = .001;
    float dt_fb = .001;
    int cmode = 1;
    float groundH = -0.404;   
    USRCMD *usrcmd = nullptr;
};

struct  TIME_PER_ITERATION
{
    int DDP_iter = 0;
    int n_bws = 0; // number of backward sweep iterations
    double time_bws =0; // backward sweep time
    int n_fit = 0; // number of iterations in line search
    double time_fit = 0; // line search time
    double time_partial = 0; // time for computing derivatives (dynamics and cost)

};

using duration_ms = std::chrono::duration<float, std::chrono::milliseconds::period>;

#ifdef TIME_BENCHMARK
extern vector<TIME_PER_ITERATION> time_ddp;
#endif //TIME_BENCHMARK

#endif // MHPC_COMPOUNDTYPES_H