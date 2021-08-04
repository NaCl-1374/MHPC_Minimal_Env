#include "SinglePhaseAbstract.h"
template <typename T>
SinglePhaseAbstract<T>::SinglePhaseAbstract(RobotBase<T> *model,
                                            CostAbstract<T> *cost,
                                            Constraint<T> *constraint,
                                            HSDDP_OPTION<T> *option,
                                            T dt, size_t N_TIMESTEPS,
                                            size_t xsize, size_t usize, size_t ysize,
                                            size_t modeidx, size_t phaseidx)
{    
    _xsize = xsize;
    _usize = usize;
    _ysize = ysize;
    set_models(model, cost, constraint, option, dt);
    set_phase_config(modeidx, phaseidx, N_TIMESTEPS);

    /* Initialize variabels */
    _V = 0;
    _dV = 0;
    _dVnext = 0;
    _Gnext.setZero(xsize);
    _Hnext.setZero(xsize,xsize);
    Ixx = DMat<T>::Identity(xsize, xsize);
    Iuu = DMat<T>::Identity(usize, usize);
    _x0.setZero(xsize);
}

template <typename T>
SinglePhaseAbstract<T>::SinglePhaseAbstract(size_t xsize, size_t usize, size_t ysize)
{
    _xsize = xsize;
    _usize = usize;
    _ysize = ysize;

    _dt = 0;
    _N_TIMESTEPS = 0;
    _modeidx = 1;
    _phaseidx = 1;
    
    /* Initialize variabels */
    _V = 0;
    _dV = 0;
    _dVnext = 0;
    _Gnext.setZero(xsize);
    _Hnext.setZero(xsize, xsize);
    Ixx = DMat<T>::Identity(xsize, xsize);
    Iuu = DMat<T>::Identity(usize, usize);
    _x0.setZero(xsize);
}

template <typename T>
void SinglePhaseAbstract<T>::set_models(RobotBase<T> *model,
                                        CostAbstract<T> *cost,
                                        Constraint<T> *constraint,
                                        HSDDP_OPTION<T> *option,
                                        T dt)
{
    _model = model;
    _cost = cost;
    _constraint = constraint;
    _option = option;
    _dt = dt;   
}

template <typename T>
void SinglePhaseAbstract<T>::set_CTG_from_nextPhase(T &dVnext, DVec<T> &Gnext, DMat<T> &Hnext)
{
    _dVnext = dVnext;
    _Gnext = Gnext;
    _Hnext = Hnext;
}

template class SinglePhaseAbstract<double>;
