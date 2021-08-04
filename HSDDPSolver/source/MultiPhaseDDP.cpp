#include "MultiPhaseDDP.h"
#include <algorithm>
#include "MHPC_CompoundTypes.h"

#ifdef TIME_BENCHMARK
#include <chrono>
using namespace std::chrono;
#endif //TIME_BENCHMARK

template <typename T>
MultiPhaseDDP<T>::MultiPhaseDDP(int n_phases,
                                HSDDP_OPTION<T> &option)
{
    _n_phases = n_phases;
    _option = option;
    _actual_cost = 0;
    _exp_cost_change = 0;
    _phases.reserve(n_phases);
    _phases.insert(_phases.begin(), n_phases, nullptr);
}

/*
  @brief: perform forward sweep for the multi-phase problem
  @params:
          eps: line search parameter (0, 1)
*/
template <typename T>
void MultiPhaseDDP<T>::forward_sweep(T eps)
{
    _actual_cost = 0;
    DVec<T> x0 = _x0;       // initial condition of first phase equal to the inital condition of the multi-phase problem
    _tconstr_violation = 0; // initialize terminal constraint violation

    for (int pidx = 0; pidx < _n_phases; pidx++)
    {
        _phases[pidx]->set_initial_condition(x0);
        _phases[pidx]->forward_sweep(eps);

        // If next state is within scope, update initial condition for the next phase
        if (pidx + 1 < _n_phases)
        {
            x0 = phase_transition(_phases[pidx], _phases[pidx + 1]);
        }
        _actual_cost += _phases[pidx]->_V;                            // add accumulated cost of each phase to total cost
        _tconstr_violation += _phases[pidx]->get_tconstr_violation(); // add up terminal constraint violation of each phase
    }
    _tconstr_violation = sqrt(_tconstr_violation);
}

/*
  @brief: perform forward sweep of dynamics propagation only for the multi-phase problem
  @params:
          eps: line search parameter (0, 1)
*/
template <typename T>
void MultiPhaseDDP<T>::forward_sweep_dynamics_only(T eps)
{
    _actual_cost = 0;
    DVec<T> x0 = _x0;       // initial condition of first phase equal to the inital condition of the multi-phase problem
    _tconstr_violation = 0; // initialize terminal constraint violation

    for (int pidx = 0; pidx < _n_phases; pidx++)
    {
        _phases[pidx]->set_initial_condition(x0);
        _phases[pidx]->forward_sweep_dynamics_only(eps);

        // If next state is within scope, update initial condition for the next phase
        if (pidx + 1 < _n_phases)
        {
            x0 = phase_transition(_phases[pidx], _phases[pidx + 1]);
        }
        _actual_cost += _phases[pidx]->_V;                            // add accumulated cost of each phase to total cost
        _tconstr_violation += _phases[pidx]->get_tconstr_violation(); // add up terminal constraint violation of each phase
    }
    _tconstr_violation = sqrt(_tconstr_violation);
}

/*
  @brief: compute dynamcis partials for the multi-phase problem
  @params:
          eps: line search parameter (0, 1)
*/
template <typename T>
void MultiPhaseDDP<T>::forward_sweep_partials_only()
{

    for (int pidx = 0; pidx < _n_phases; pidx++)
    {
        _phases[pidx]->forward_sweep_partials_only();
    }
}

/*
  @brief: perform backward sweep for the multi-phase problem
  @params: 
          regularization: regularization parameter
  @return: success of backward sweep
*/
template <typename T>
bool MultiPhaseDDP<T>::backward_sweep(T regularization)
{
    bool success = true;
    dVnext = 0;
    Gnext.setZero(_phases[_n_phases - 1]->_xsize);
    Hnext.setZero(_phases[_n_phases - 1]->_xsize, _phases[_n_phases - 1]->_xsize);

    for (int pidx = _n_phases - 1; pidx >= 0; pidx--)
    {
        // If not the last phase, udpate CTG information back propagated from next phase
        if (pidx + 1 < _n_phases)
        {
            impact_aware_step(_phases[pidx], _phases[pidx + 1], dVnext, Gnext, Hnext);
        }
        _phases[pidx]->set_CTG_from_nextPhase(dVnext, Gnext, Hnext);

        // If backward sweep of current phase fails, break and return false
        if (!_phases[pidx]->backward_sweep(regularization))
        {
            success = false;
            break;
        }
    }

    if (success)
        _exp_cost_change = _phases[0]->_dV;
    return success;
}

template <typename T>
int MultiPhaseDDP<T>::forward_iteration()
{
    T eps = 1;
    T cost_prev = _actual_cost;
    int iter = 1;
    while (eps > pow(0.1, 10))
    {
        forward_sweep_dynamics_only(eps);
#ifdef DEBUG
        printf("\t eps=%.3e \t actual change in cost=%.3e \t expeced change in cost=%.3e\n",
               eps, _actual_cost - cost_prev, _option.gamma * eps * (1 - eps / 2) * _exp_cost_change);
#endif
        if (_actual_cost <= cost_prev + _option.gamma * eps * (1 - eps / 2) * _exp_cost_change)
        {
            break;
        }

        eps *= _option.alpha;
        iter++;
    }
    return iter;
}

template <typename T>
void MultiPhaseDDP<T>::solve()
{
    int iter_AL = 1;
    int iter_DDP_abs = 1;
    T cost_prev;
    bool bws_success;
    T update_penalty = _option.update_penalty;
    bool ReB_active = _option.ReB_active;

#ifdef TIME_BENCHMARK
    time_ddp.clear();
    double time_partial = 0;
    TIME_PER_ITERATION time_per_iter;
    auto start = high_resolution_clock::now();
    auto stop = high_resolution_clock::now();
    auto duration = duration_ms(stop - start);
#endif

    while (iter_AL <= _option.max_AL_iter)
    {
#ifdef DEBUG
        printf("outer loop iteration %d \n", iter_AL);
#endif

        _option.ReB_active = ReB_active;
        if ((_tconstr_violation > 0.05) || 1 == iter_AL)
        /* If terminal constraint violation is too large, turn path constraint off */
        {
            _option.ReB_active = 0;
        }
#ifdef TIME_BENCHMARK
        start = high_resolution_clock::now();
#endif
        forward_sweep(0);
#ifdef TIME_BENCHMARK
        stop = high_resolution_clock::now();
        duration = duration_ms(stop - start);
        time_partial = duration.count();
#endif
        update_nominal_trajectory();
        int iter_DDP = 1;
        T regularization = 0;
        while (iter_DDP <= _option.max_DDP_iter)
        {
#ifdef DEBUG
            printf("\t inner loop iteration %d \n", iter_DDP);
#endif
            cost_prev = _actual_cost;

            bws_success = false;
            int bws_iter = 1;

#ifdef TIME_BENCHMARK
            auto start = high_resolution_clock::now();
#endif
            while (!bws_success)
            {
#ifdef DEBUG
                printf("\t regularization = %f\n", regularization);
#endif
                bws_success = backward_sweep(regularization);
                if (bws_success)
                    break;

                regularization = std::max(regularization * _option.update_regularization, T (1e-03));
                bws_iter++;

                if (regularization > 1000)
                {
                    printf("Regularization term exceeds maximum value! \n");
                    printf("Optimization terminates! \n");
                    return;
                }
            }
#ifdef TIME_BENCHMARK
            stop = high_resolution_clock::now();
            duration = duration_ms(stop - start);
            time_per_iter.n_bws = bws_iter;
            time_per_iter.time_bws = duration.count();
            time_per_iter.time_partial = time_partial;
            time_per_iter.DDP_iter = iter_DDP_abs;
#endif

            regularization = regularization / 20;
            if (regularization < 1e-06)
            {
                regularization = 0;
            }
#ifdef TIME_BENCHMARK
            start = high_resolution_clock::now();
#endif
            int fit_iter = forward_iteration(); // Performe line search
#ifdef TIME_BENCHMARK
            stop = high_resolution_clock::now();
            duration = duration_ms(stop - start);
            time_per_iter.n_fit = fit_iter;
            time_per_iter.time_fit = duration.count();
            time_ddp.push_back(time_per_iter);

#endif

            update_nominal_trajectory();
            if (cost_prev - _actual_cost < _option.DDP_thresh)
                break; // If DDP converges, break innter loop

#ifdef TIME_BENCHMARK
            start = high_resolution_clock::now();
#endif
            forward_sweep_partials_only();
#ifdef TIME_BENCHMARK
            stop = high_resolution_clock::now();
            duration = duration_ms(stop - start);
            time_partial = duration.count();
#endif

            iter_DDP++;
            iter_DDP_abs++;
        }

        _option.update_penalty = update_penalty;
        if (_tconstr_violation < 0.03)
        {
            _option.update_penalty = 0;
        }
        update_AL_ReB_param();

#ifdef DEBUG
        printf("terminal constraint violation = %.3e \n", _tconstr_violation);
#endif

        if (_tconstr_violation < _option.AL_thresh)
            break;

        iter_AL++;
    }
}

/*
  @brief: perform impact-aware backward DDP step
  @params: 
          cPhase: current phase
          nPhase: next phase
          dV, G, H: expected cost reduction, gradient and hessian of value function propaged from next phase through phase transition
*/

template <typename T>
void MultiPhaseDDP<T>::impact_aware_step(SinglePhaseAbstract<T> *cPhase, SinglePhaseAbstract<T> *nPhase, T &dV, DVec<T> &G, DMat<T> &H)
{
    DVec<T> Gprime(nPhase->_xsize);
    DMat<T> Hprime(nPhase->_xsize, nPhase->_xsize);

    Gprime.setZero();
    Hprime.setZero();

    // Get CTG information from next phase
    nPhase->get_CTG_currentPhase(dV, Gprime, Hprime);

    // If current phase is full model, perform impact-aware DDP step
    if (ModelType::FULL == cPhase->get_model_type())
    {
        MatMN<T, xsize_WB, xsize_WB> Px;
        VecM<T, xsize_WB> xend;
        Px.setZero();

        xend = cPhase->get_terminal_state();

        cPhase->_model->resetmap_par(xend, Px, cPhase->get_modeidx());

        switch (nPhase->get_model_type())
        {
        case ModelType::FULL:
            G = Px.transpose() * Gprime;
            H = Px.transpose() * Hprime * Px;
            break;

        case ModelType::RBD: // If next phase is simple model, include the state projection into impact-aware step
            G = Px.transpose() * _stateProj.transpose() * Gprime;
            H = Px.transpose() * _stateProj.transpose() * Hprime * _stateProj * Px;
            break;
        }
    }
    // If current phase is simple model, the transition is smooth
    else if (ModelType::RBD == cPhase->get_model_type())
    {
        G = Gprime;
        H = Hprime;
    }
}

/*
  @brief: propogates state from the end of current phase to the begining of next phase
  @params: 
          cPhase: current phase
          nPhase: next phase
  @return: initial condition of next phase
*/
template <typename T>
DVec<T> MultiPhaseDDP<T>::phase_transition(SinglePhaseAbstract<T> *cPhase, SinglePhaseAbstract<T> *nPhase)
{

    // If current phase is full model, use resetmap
    if (ModelType::FULL == cPhase->get_model_type())
    {
        VecM<T, xsize_WB> xend_cphase, x_interm;
        VecM<T, ysize_WB> y;
        xend_cphase = cPhase->get_terminal_state();
        x_interm.setZero();
        y.setZero();

        cPhase->_model->resetmap(xend_cphase, x_interm, y, cPhase->get_modeidx());

        switch (nPhase->get_model_type())
        {
        case ModelType::FULL: // If next phase is full model, pass x_interm to x0_nphase
            return x_interm;

        case ModelType::RBD: // If next phase is simple model, perform state projection
            return _stateProj * x_interm;
        }
    }
    // If current phase is RBD model, use identity map
    else
    {
        return cPhase->get_terminal_state();
    }
}

template <typename T>
void MultiPhaseDDP<T>::empty_bag()
{
    for (int pidx = 0; pidx < _n_phases; pidx++)
    {
        if (!_phases[pidx])
        {
            _phases[pidx]->empty_bag();
        }
    }
}
template <typename T>
void MultiPhaseDDP<T>::update_nominal_trajectory()
{
    for (int pidx = 0; pidx < _n_phases; pidx++)
    {
        if (nullptr != _phases[pidx])
        {
            _phases[pidx]->update_nominal_trajectory();
        }
    }
}

template <typename T>
void MultiPhaseDDP<T>::update_AL_ReB_param()
{
    for (int pidx = 0; pidx < _n_phases; pidx++)
    {
        if (nullptr != _phases[pidx])
        {
            _phases[pidx]->udpate_AL_ReB_Param();
        }
    }
}

/*
  @brief: warm start the multi-phase trajecotry.
          full-model phase is warm started with heuristic controller
          simple-model phase is initialized with zeros
*/
template <typename T>
void MultiPhaseDDP<T>::warmstart()
{
    for (int pidx = 0; pidx < _n_phases; pidx++)
    {
        if (nullptr != _phases[pidx])
        {
            _phases[pidx]->set_data_to_zero();
        }
    }
}

template class MultiPhaseDDP<double>;
