#include"PlanarFloatingBase.h"
#include"FBDynamics.h"
#include"FBDynamics_par.h"
#include"CasadiGen.h"

template<typename T>
void PlanarFloatingBase<T>::dynamics(VecM<T,xsize_FB> &x, VecM<T,usize_FB> &u, VecM<T, xsize_FB> &x_next, VecM<T,ysize_FB>&y, int mode)
{
    switch (mode)
    {
    case 1:
        _contact_state << 0, 1;
        break;
    case 2:
        _contact_state << 0, 0;
        break;
    case 3:
        _contact_state << 1, 0;
        break;
    case 4:
        _contact_state << 0, 0;
        break;                
    }
    _xdot.setZero();
    vector<T *> arg = {x.data(), u.data(), _foothold.data(), (T *)_contact_state.data()};
    vector<T *> res = {_xdot.data()};

    casadi_interface(arg, res, x.size(), FBDynamics, FBDynamics_sparsity_out, FBDynamics_work);

    y.setZero();
    x_next = x + _xdot* this->_dt;
}

template<typename T>
void PlanarFloatingBase<T>::resetmap(VecM<T,xsize_FB> &x, VecM<T,xsize_FB> &x_next, VecM<T,ysize_FB>&y, int mode)
// void PlanarFloatingBase<T>::resetmap(DVec<T> &x, DVec<T> &x_next, DVec<T> &y, int mode = 1)

{
    x_next = x;
    y.setZero();
}

template<typename T>
void PlanarFloatingBase<T>::dynamics_par(VecM<T,xsize_FB>& x, VecM<T,usize_FB> &u, MatMN<T,xsize_FB,xsize_FB> &A, MatMN<T,xsize_FB,usize_FB> &B, 
                              MatMN<T,ysize_FB, xsize_FB> &C, MatMN<T, ysize_FB, usize_FB> &D, int mode)
{
    switch (mode)
    {
    case 1:
        _contact_state << 0, 1;
        break;
    case 2:
        _contact_state << 0, 0;
        break;
    case 3:
        _contact_state << 1, 0;
        break;
    case 4:
        _contact_state << 0, 0;
        break;                
    }
    _Ac.setZero();
    _Bc.setZero();
    vector<T *> arg = {x.data(), u.data(), _foothold.data(), (T *)_contact_state.data()};
    vector<T *> res = {_Ac.data(), _Bc.data()};

    casadi_interface(arg, res, _Ac.size(), FBDynamics_par, FBDynamics_par_sparsity_out, FBDynamics_par_work);
    A = MatMN<T,xsize_FB,xsize_FB>::Identity() + _Ac* this->_dt;
    B = _Bc * this->_dt;
    C.setZero();
    D.setZero();
}

template<typename T>
void PlanarFloatingBase<T>::resetmap_par(VecM<T, xsize_FB> &x, MatMN<T,xsize_FB, xsize_FB> &Px, int mode)
{
    Px.setIdentity();
}

template<typename T>
void PlanarFloatingBase<T>::plan_foothold(DVec<T> &x,T stance_time, int mode) 
{
    _foothold= _foothold_planner->get_foothold_location(x.head(3), stance_time, mode);
}

template class PlanarFloatingBase<double>;