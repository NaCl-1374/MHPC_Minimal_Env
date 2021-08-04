#ifndef MULTIPHASEDDP_H
#define MULTIPHASEDDP_H

#include <vector>
#include "SinglePhase.h"
#include "PlanarRobot.h"

template <typename T>
class MultiPhaseDDP
{
// protected:
public:
    // SinglePhaseAbstract<T> **_phases = nullptr;
    std::vector<SinglePhaseAbstract<T> *> _phases;

public:
    MultiPhaseDDP() {} // default constructor
    MultiPhaseDDP(int n_phases,
                  HSDDP_OPTION<T> &option);        
    virtual ~MultiPhaseDDP() = default;
    
    void set_option(HSDDP_OPTION<T> &option) {_option = option;}
    
    template <typename Derived>
    void set_initial_condition(const Eigen::MatrixBase<Derived> &x0) { _x0 = x0; }

    virtual void warmstart();

    void forward_sweep(T eps);

    void forward_sweep_dynamics_only(T eps);

    void forward_sweep_partials_only();

    bool backward_sweep(T regularization);

    int forward_iteration();

    void solve();

    void impact_aware_step(SinglePhaseAbstract<T> *cPhase, SinglePhaseAbstract<T> *nPhase, T &, DVec<T> &, DMat<T> &);

    DVec<T> phase_transition(SinglePhaseAbstract<T> *cPhase, SinglePhaseAbstract<T> *nPhase);

    virtual void update_nominal_trajectory();

    void empty_bag();

    void initialize_AL_ReB_param();

    void update_AL_ReB_param();

public:
    int _n_phases;
    DMat<T> _stateProj; // low-rank state projection from full model to simple model

    T _actual_cost;
    T _exp_cost_change;

    T _tconstr_violation;

    HSDDP_OPTION<T> _option;
    DVec<T> _x0;

private:
    DVec<T> Gnext;
    DMat<T> Hnext;
    T dVnext;
};

#endif // MULTIPHASEDDP_H