
#include <assert.h>
#include <math.h>
#include "SinglePhase.h"
#include <iostream>

template <typename T, size_t XSIZE, size_t USIZE, size_t YSIZE>
SinglePhase<T, XSIZE, USIZE, YSIZE>::SinglePhase(RobotBase<T> *model,
                                                 CostAbstract<T> *cost,
                                                 Constraint<T> *constraint,
                                                 HSDDP_OPTION<T> *option,
                                                 T dt, size_t N_TIMESTEPS,
                                                 size_t modeidx, size_t phaseidx)
    : SinglePhaseAbstract<T>(model,
                             cost,
                             constraint,
                             option,
                             dt, N_TIMESTEPS,
                             XSIZE, USIZE, YSIZE,
                             modeidx, phaseidx)
{
    initialization();
}

template <typename T, size_t XSIZE, size_t USIZE, size_t YSIZE>
void SinglePhase<T, XSIZE, USIZE, YSIZE>::initialization()
{
    _num_pconstr = this->_constraint->get_num_pconstraint(this->_modeidx);
    _num_tconstr = this->_constraint->get_num_tconstraint(this->_modeidx);   
    
    if (_pconstr == nullptr)
        delete [] _pconstr;
    if (_tconstr == nullptr)
        delete [] _tconstr;
    _pconstr = new IneqConstrStruct<T, XSIZE, USIZE, YSIZE>[_num_pconstr];
    _tconstr = new TConstrStruct<T, XSIZE>[_num_tconstr];
    B.setZero(_num_pconstr);
    Bz.setZero(_num_pconstr);
    Bzz.setZero(_num_pconstr);
    // Initialize AL and ReB parameter
    initialize_AL_ReB_Params();
}

template <typename T, size_t XSIZE, size_t USIZE, size_t YSIZE>
void SinglePhase<T, XSIZE, USIZE, YSIZE>::set_data(ModelState<T, XSIZE, USIZE, YSIZE> *ms_act, ModelState<T, XSIZE, USIZE, YSIZE> *ms_nom,
                                                   DynDerivative<T, XSIZE, USIZE, YSIZE> *dynpar,
                                                   RCostStruct<T, XSIZE, USIZE, YSIZE> *rcost,
                                                   CostToGoStruct<T, XSIZE, USIZE> *CTG,
                                                   TCostStruct<T, XSIZE> *tcost,
                                                   ModelState<T, XSIZE, USIZE, YSIZE> *ms_ref)
{
    _ms_act = ms_act;
    _ms_nom = ms_nom;
    _dynpar = dynpar;
    _rcost = rcost;
    _CTG = CTG;
    _tcost = tcost;
    _ms_ref = ms_ref;
}

template <typename T, size_t XSIZE, size_t USIZE, size_t YSIZE>
void SinglePhase<T, XSIZE, USIZE, YSIZE>::forward_sweep(T eps)
{
    this->_V = 0;
    int k;
    _ms_act[0].x = this->_x0;
    // If rigid-body model, plan the foothold location first
    if (ModelType::RBD == this->_model->get_model_type())
    {
        this->_model->plan_foothold(this->_x0, this->_dt * this->_N_TIMESTEPS, this->_modeidx);
    }

    for (k = 0; k < this->_N_TIMESTEPS - 1; k++)
    {
        /* update control */
        _ms_act[k].u = _ms_nom[k].u + eps * _CTG[k].du + _CTG[k].K * (_ms_act[k].x - _ms_nom[k].x);

        /* run dynamics */
        this->_model->dynamics(_ms_act[k].x, _ms_act[k].u, _ms_act[k + 1].x, _ms_act[k].y, this->_modeidx);

        /* compuate dynamics derivatives */
        this->_model->dynamics_par(_ms_act[k].x, _ms_act[k].u, _dynpar[k].A, _dynpar[k].B, _dynpar[k].C, _dynpar[k].D, this->_modeidx);

        /* compute running cost */
        this->_cost->running_cost(_ms_act[k], _ms_ref[k], _rcost[k], this->_modeidx);

        /* compute running cost derivatives */
        this->_cost->running_cost_par(_ms_act[k], _ms_ref[k], _rcost[k], this->_modeidx);

        /* compute inequality constraint */
        this->_constraint->path_constraint(_ms_act[k], _pconstr, this->_modeidx);

        /* update running cost with ineq constraint encoded by ReB functions */
        if ((this->_option->ReB_active) && (!_param.reb_empty))
            update_running_cost_with_pconstr(_rcost[k], _pconstr, _param.delta, _param.eps_ReB, CALC_DYN_AND_PAR);

        this->_V += _rcost[k].l;
    }

    /* compute terminal cost */
    this->_cost->terminal_cost(_ms_act[this->_N_TIMESTEPS - 1], _ms_ref[this->_N_TIMESTEPS - 1], *_tcost, this->_modeidx);

    /* compute terminal cost partial information */
    this->_cost->terminal_cost_par(_ms_act[this->_N_TIMESTEPS - 1], _ms_ref[this->_N_TIMESTEPS - 1], *_tcost, this->_modeidx);

    /* compute terminal constraint */
    this->_constraint->terminal_constraint(_ms_act[this->_N_TIMESTEPS - 1], _tconstr, this->_modeidx);

    /* update terminal cost with terminal constraint encoded by AL functions */
    if ((this->_option->AL_active) & (!_param.al_empty))
        update_terminal_cost_with_tconstr(*_tcost, _tconstr, _param.sigma, _param.lambda, CALC_DYN_AND_PAR);

    this->_V += _tcost->Phi;
}

template <typename T, size_t XSIZE, size_t USIZE, size_t YSIZE>
void SinglePhase<T, XSIZE, USIZE, YSIZE>::forward_sweep_dynamics_only(T eps)
{
    this->_V = 0;
    int k;
    _ms_act[0].x = this->_x0;
    // If rigid-body model, plan the foothold location first
    if (ModelType::RBD == this->_model->get_model_type())
    {
        this->_model->plan_foothold(this->_x0, this->_dt * this->_N_TIMESTEPS, this->_modeidx);
    }

    for (k = 0; k < this->_N_TIMESTEPS - 1; k++)
    {
        _ms_act[k].u = _ms_nom[k].u + eps * _CTG[k].du + _CTG[k].K * (_ms_act[k].x - _ms_nom[k].x);         //update control
        this->_model->dynamics(_ms_act[k].x, _ms_act[k].u, _ms_act[k + 1].x, _ms_act[k].y, this->_modeidx); //run dynamics
        this->_cost->running_cost(_ms_act[k], _ms_ref[k], _rcost[k], this->_modeidx);                       // compute running cost
        this->_constraint->path_constraint(_ms_act[k], _pconstr, this->_modeidx);                           // compute inequality constraint
        if ((this->_option->ReB_active) && (!_param.reb_empty))                                             // update running cost with ineq constraint
            update_running_cost_with_pconstr(_rcost[k], _pconstr, _param.delta, _param.eps_ReB, CALC_DYNAMICS_ONLY);
        this->_V += _rcost[k].l;
    }

    this->_cost->terminal_cost(_ms_act[this->_N_TIMESTEPS - 1], _ms_ref[this->_N_TIMESTEPS - 1], *_tcost, this->_modeidx); // compute terminal cost
    this->_constraint->terminal_constraint(_ms_act[this->_N_TIMESTEPS - 1], _tconstr, this->_modeidx);                     // compute terminal constraint
    if ((this->_option->AL_active) & (!_param.al_empty))                                                                   // update terminal cost with terminal constraint encoded by AL functions
        update_terminal_cost_with_tconstr(*_tcost, _tconstr, _param.sigma, _param.lambda, CALC_DYNAMICS_ONLY);
    this->_V += _tcost->Phi;
}

template <typename T, size_t XSIZE, size_t USIZE, size_t YSIZE>
void SinglePhase<T, XSIZE, USIZE, YSIZE>::forward_sweep_partials_only()
{
    // If rigid-body model, plan the foothold location first
    if (ModelType::RBD == this->_model->get_model_type())
    {
        this->_model->plan_foothold(this->_x0, this->_dt * this->_N_TIMESTEPS, this->_modeidx);
    }

    for (int k = 0; k < this->_N_TIMESTEPS - 1; k++)
    {
        /* compuate dynamics derivatives */
        this->_model->dynamics_par(_ms_act[k].x, _ms_act[k].u, _dynpar[k].A, _dynpar[k].B, _dynpar[k].C, _dynpar[k].D, this->_modeidx);

        /* compute running cost derivatives */
        this->_cost->running_cost_par(_ms_act[k], _ms_ref[k], _rcost[k], this->_modeidx);

        /* compute inequality constraint */
        this->_constraint->path_constraint(_ms_act[k], _pconstr, this->_modeidx);

        /* update running cost with ineq constraint encoded by ReB functions */
        if ((this->_option->ReB_active) && (!_param.reb_empty))
            update_running_cost_with_pconstr(_rcost[k], _pconstr, _param.delta, _param.eps_ReB, CALC_PARTIALS_ONLY);
    }

    /* compute terminal cost partial information */
    this->_cost->terminal_cost_par(_ms_act[this->_N_TIMESTEPS - 1], _ms_ref[this->_N_TIMESTEPS - 1], *_tcost, this->_modeidx);

    /* compute terminal constraint */
    this->_constraint->terminal_constraint(_ms_act[this->_N_TIMESTEPS - 1], _tconstr, this->_modeidx);

    /* update terminal cost with terminal constraint encoded by AL functions */
    if ((this->_option->AL_active) & (!_param.al_empty))
        update_terminal_cost_with_tconstr(*_tcost, _tconstr, _param.sigma, _param.lambda, CALC_PARTIALS_ONLY);
}

template <typename T, size_t XSIZE, size_t USIZE, size_t YSIZE>
bool SinglePhase<T, XSIZE, USIZE, YSIZE>::backward_sweep(T regularization)
{
    bool success = true;

    /*    _Gnext, _Hnext, _dVnext
    */
    this->_dV = this->_dVnext;
    _CTG[this->_N_TIMESTEPS - 1].G = _tcost->Phix + this->_Gnext;
    _CTG[this->_N_TIMESTEPS - 1].H = _tcost->Phixx + this->_Hnext;

    for (int k = this->_N_TIMESTEPS - 2; k >= 0; k--)
    {
        // compute Q function approximation
        _CTG[k].compute_Qfunction(_rcost[k], _dynpar[k], _CTG[k + 1].G, _CTG[k + 1].H);

        // regularizatoin
        _CTG[k].Qxx += (this->Ixx * regularization);
        _CTG[k].Quu += (this->Iuu * regularization);

        this->Quu_chol = this->Quu_chol.compute(_CTG[k].Quu - this->Iuu * pow(0.1, 9)); // Cholesky decomposition of Quu

        // If Quu not PSD, break and return false
        if (!this->Quu_chol.isPositive())
        {
            success = false;
            break;
        }

        // compute value function approximation and sum up expected cost change at each time step
        this->_dV += _CTG[k].valuefunction_update();
    }

    return success;
}

template <typename T, size_t XSIZE, size_t USIZE, size_t YSIZE>
void SinglePhase<T, XSIZE, USIZE, YSIZE>::update_running_cost_with_pconstr(
    RCostStruct<T, XSIZE, USIZE, YSIZE> &Rcost,
    IneqConstrStruct<T, XSIZE, USIZE, YSIZE>* Ineq,
    DVec<T> &delta, DVec<T> &eps, int calcflag)
{
    if (Ineq == nullptr)
        return;

    B.setZero();
    Bz.setZero();
    Bzz.setZero();
    reduced_barrier(Ineq, delta, B, Bz, Bzz);

    for (size_t idx = 0; idx < B.size(); idx++)
    {
        if (calcflag == CALC_DYNAMICS_ONLY || calcflag == CALC_DYN_AND_PAR)
        {
            Rcost.l += eps(idx) * B(idx) * this->_dt;
        }

        if (calcflag == CALC_PARTIALS_ONLY || calcflag == CALC_DYN_AND_PAR)
        {
            Rcost.lx += eps(idx) * Bz(idx) * Ineq[idx].gx * this->_dt;
            Rcost.lu += eps(idx) * Bz(idx) * Ineq[idx].gu * this->_dt;
            Rcost.ly += eps(idx) * Bz(idx) * Ineq[idx].gy * this->_dt;
            Rcost.lxx += eps(idx) * (Ineq[idx].gx * Bzz(idx) * Ineq[idx].gx.transpose() + Bz(idx) * Ineq[idx].gxx) * this->_dt;
            Rcost.luu += eps(idx) * (Ineq[idx].gu * Bzz(idx) * Ineq[idx].gu.transpose() + Bz(idx) * Ineq[idx].guu) * this->_dt;
            Rcost.lyy += eps(idx) * (Ineq[idx].gy * Bzz(idx) * Ineq[idx].gy.transpose() + Bz(idx) * Ineq[idx].gyy) * this->_dt;
        }
    }
}

template <typename T, size_t XSIZE, size_t USIZE, size_t YSIZE>
void SinglePhase<T, XSIZE, USIZE, YSIZE>::update_running_cost_with_smooth()
{
}

template <typename T, size_t XSIZE, size_t USIZE, size_t YSIZE>
void SinglePhase<T, XSIZE, USIZE, YSIZE>::update_terminal_cost_with_tconstr(
    TCostStruct<T, XSIZE> &Tcost, TConstrStruct<T, XSIZE>* Tconstr, T sigma, DVec<T> &lambda, int calcflag)
{
    if (Tconstr == nullptr)
        return;

    for (size_t idx = 0; idx < _num_tconstr; idx++)
    {
        if (calcflag == CALC_DYNAMICS_ONLY || calcflag == CALC_DYN_AND_PAR)
        {
            Tcost.Phi += 50 * (pow(sigma * Tconstr[idx].h / 2, 2) + lambda[idx] * Tconstr[idx].h);
        }
        if (calcflag == CALC_DYNAMICS_ONLY || calcflag == CALC_DYN_AND_PAR)
        {
            Tcost.Phix += 50 * (sigma * sigma / 2 * Tconstr[idx].hx * Tconstr[idx].h + lambda[idx] * Tconstr[idx].hx);
            Tcost.Phixx += 50 * (sigma * sigma / 2 * (Tconstr[idx].hx * Tconstr[idx].hx.transpose() + Tconstr[idx].h * Tconstr[idx].hxx) + lambda[idx] * Tconstr[idx].hxx);
        }
    }
}

template <typename T, size_t XSIZE, size_t USIZE, size_t YSIZE>
void SinglePhase<T, XSIZE, USIZE, YSIZE>::empty_bag()
{
    if (nullptr != _ms_act)
        std::memset(_ms_act, 0, this->_N_TIMESTEPS * sizeof(ModelState<T, XSIZE, USIZE, YSIZE>));
    if (nullptr != _CTG)
        std::memset(_CTG, 0, this->_N_TIMESTEPS * sizeof(CostToGoStruct<T, XSIZE, USIZE>));
}

template <typename T, size_t XSIZE, size_t USIZE, size_t YSIZE>
void SinglePhase<T, XSIZE, USIZE, YSIZE>::set_data_to_zero()
{
    if (nullptr != _ms_act)
        std::memset(_ms_act, 0, this->_N_TIMESTEPS * sizeof(ModelState<T, XSIZE, USIZE, YSIZE>));
    if (nullptr != _ms_nom)
        std::memset(_ms_act, 0, this->_N_TIMESTEPS * sizeof(ModelState<T, XSIZE, USIZE, YSIZE>));
    if (nullptr != _CTG)
        std::memset(_CTG, 0, this->_N_TIMESTEPS * sizeof(CostToGoStruct<T, XSIZE, USIZE>));
}

template <typename T, size_t XSIZE, size_t USIZE, size_t YSIZE>
void SinglePhase<T, XSIZE, USIZE, YSIZE>::reduced_barrier(
    IneqConstrStruct<T, XSIZE, USIZE, YSIZE> *Ineq, DVec<T> &delta, DVec<T> &b, DVec<T> &bz, DVec<T> &bzz)
{
    int k = 2; // order of approximating polynomial
    for (size_t idx = 0; idx < b.size(); idx++)
    {
        if (Ineq[idx].g > delta[idx])
        {
            b[idx] = -log(Ineq[idx].g);
            bz[idx] = -1.0 / Ineq[idx].g;
            bzz[idx] = pow(Ineq[idx].g, -2);
        }
        else
        {
            b[idx] = (double)(k - 1) / k * (pow((Ineq[idx].g - k * delta[idx]) / ((k - 1) * delta[idx]), k) - 1) - log(delta[idx]);
            bz[idx] = pow((Ineq[idx].g - k * delta[idx]) / ((k - 1) * delta[idx]), (k - 1)) / delta[idx];
            bzz[idx] = pow((Ineq[idx].g - k * delta[idx]) / ((k - 1) * delta[idx]), k - 2);
        }
    }
}

template <typename T, size_t XSIZE, size_t USIZE, size_t YSIZE>
void SinglePhase<T, XSIZE, USIZE, YSIZE>::get_CTG_currentPhase(T &dV, DVec<T> &G, DMat<T> &H)
{
    dV = this->_dV;
    G = _CTG[0].G;
    H = _CTG[0].H;
}

template <typename T, size_t XSIZE, size_t USIZE, size_t YSIZE>
void SinglePhase<T, XSIZE, USIZE, YSIZE>::initialize_AL_ReB_Params()
{
    this->_constraint->get_AL_REB_PARAMS(_param, this->_modeidx);
}

template <typename T, size_t XSIZE, size_t USIZE, size_t YSIZE>
void SinglePhase<T, XSIZE, USIZE, YSIZE>::udpate_AL_ReB_Param()
{
    for (size_t eqidx = 0; eqidx < _param.lambda.size(); eqidx++)
    {
        _param.lambda[eqidx] += _param.sigma * _tconstr[eqidx].h;
    }
    _param.sigma *= this->_option->update_penalty;

    if (this->_option->ReB_active)
    {
        for (size_t ineqidx = 0; ineqidx < _param.delta.size(); ineqidx++)
        {
            _param.delta[ineqidx] *= this->_option->update_relax;
            if (_param.delta[ineqidx] < _param.delta_min[ineqidx])
            {
                _param.delta[ineqidx] = _param.delta_min[ineqidx];
            }
            _param.eps_ReB[ineqidx] *= this->_option->update_ReB;
        }
    }
}


template <typename T, size_t XSIZE, size_t USIZE, size_t YSIZE>
void SinglePhase<T, XSIZE, USIZE, YSIZE>::update_nominal_trajectory()
{
    std::memcpy(_ms_nom, _ms_act, this->_N_TIMESTEPS * sizeof(ModelState<T, XSIZE, USIZE, YSIZE>));
}

template <typename T, size_t XSIZE, size_t USIZE, size_t YSIZE>
void* SinglePhase<T, XSIZE, USIZE, YSIZE>::get_nominal_ms_ptr()
{
    return _ms_nom;
}

template <typename T, size_t XSIZE, size_t USIZE, size_t YSIZE>
void* SinglePhase<T, XSIZE, USIZE, YSIZE>::get_CTG_info_ptr()
{
    return _CTG;
}


template <typename T, size_t XSIZE, size_t USIZE, size_t YSIZE>
SinglePhase<T, XSIZE, USIZE, YSIZE>::~SinglePhase()
{
}

template class SinglePhase<double, 14, 4, 4>;
template class SinglePhase<double, 6, 4, 4>;
