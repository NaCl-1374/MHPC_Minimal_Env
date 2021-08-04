#ifndef SINGLEPHASE_ABSTRACT_H
#define SINGLEPHASE_ABSTRACT_H

#include "MHPC_CPPTypes.h"
#include "MHPC_CompoundTypes.h"
#include "PlanarRobot.h"
#include "CostBase.h"
#include "ConstraintsBase.h"

template<typename> class MultiPhaseDDP; //forward declaration of MultiPhaseDDP class template

template<typename T>
class SinglePhaseAbstract
{

private:    
    friend class MultiPhaseDDP<T>; 
    SinglePhaseAbstract<T> * _nextPhase = nullptr; // ponter to next phase


public:
    SinglePhaseAbstract(RobotBase<T> *model,
                     CostAbstract<T> *cost,
                     Constraint<T> *constraint,
                     HSDDP_OPTION<T> *option,
                     T dt,  size_t N_TIMESTEPS,                     
                     size_t xsize, size_t usize, size_t ysize,
                     size_t modeidx, size_t phaseidx);

    SinglePhaseAbstract(size_t xsize, size_t usize, size_t ysize);
    virtual ~SinglePhaseAbstract() = default;

public:
    void set_models(RobotBase<T> *model,
                    CostAbstract<T> *cost,
                    Constraint<T> *constraint,
                    HSDDP_OPTION<T> *option,
                    T dt);

    virtual void set_data(ModelState<T, xsize_WB, usize_WB, ysize_WB>*, 
                          ModelState<T, xsize_WB, usize_WB, ysize_WB>*,
                          DynDerivative<T, xsize_WB, usize_WB, ysize_WB> *,
                          RCostStruct<T, xsize_WB, usize_WB, ysize_WB> *,
                          CostToGoStruct<T, xsize_WB, usize_WB> *,
                          TCostStruct<T, xsize_WB> *,
                          ModelState<T, xsize_WB, usize_WB, ysize_WB> *) {}

    virtual void set_data(ModelState<T, xsize_FB, usize_FB, ysize_FB>*, 
                          ModelState<T, xsize_FB, usize_FB, ysize_FB>*,
                          DynDerivative<T, xsize_FB, usize_FB, ysize_FB> *,
                          RCostStruct<T, xsize_FB, usize_FB, ysize_FB> *,
                          CostToGoStruct<T, xsize_FB, usize_FB> *,
                          TCostStruct<T, xsize_FB> *,
                          ModelState<T, xsize_FB, usize_FB, ysize_FB> *) {}                          

    void set_phase_config(int modeidx, int phaseidx, int N_TIMESTEPS) {_modeidx=modeidx; _phaseidx = phaseidx; _N_TIMESTEPS = N_TIMESTEPS;}

    virtual void set_reference(ModelState<T,14,4,4> *) {}

    virtual void set_reference(ModelState<T,6,4,4> *) {}

    virtual void initialization(){}

    virtual void set_data_to_zero() {}
    
    size_t get_modeidx() { return _modeidx; }

    ModelType get_model_type() {return _model->get_model_type();}

    void set_CTG_from_nextPhase(T& dVnext, DVec<T> &Gnext, DMat<T> &Hnext);

    template<typename Derived>
    void set_initial_condition(const Eigen::DenseBase<Derived> &x0) {_x0 = x0;}

    virtual void get_CTG_currentPhase(T &dV, DVec<T> &G, DMat<T> &H) {}    

    virtual DVec<T> get_terminal_state() {}

    virtual void* get_nominal_ms_ptr(){return nullptr;}

    virtual void* get_CTG_info_ptr(){return nullptr;}

    virtual void forward_sweep(T eps) = 0; // compute both dynamics propagation and dynamics partials

    virtual void forward_sweep_dynamics_only(T eps) = 0; // only compute dynamics propagation

    virtual void forward_sweep_partials_only() = 0; // only compute dynamics partials

    virtual bool backward_sweep(T regularization) = 0;
    
    virtual void update_nominal_trajectory() = 0;

    virtual void udpate_AL_ReB_Param() = 0;

    virtual void initialize_AL_ReB_Params() = 0;

    virtual void empty_bag()=0;

    virtual T get_tconstr_violation()=0;

    void set_next_phase(SinglePhaseAbstract<T> * nextPhase) {this->_nextPhase = nextPhase;}

    SinglePhaseAbstract<T> * get_next_phase() {return _nextPhase;}


public:
    // EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    RobotBase<T> *_model = nullptr;
    CostAbstract<T> *_cost = nullptr;
    Constraint<T> *_constraint = nullptr;
    HSDDP_OPTION<T> *_option =nullptr;

    int _modeidx = 1;
    int _phaseidx = 1;
    T _dt = 0;
    size_t _N_TIMESTEPS = 0;

    T _V; //actual cost for the current phase only
    T _dV; //expected cost change

    // value function information from next phase
    T _dVnext;
    DVec<T> _Gnext;
    DMat<T> _Hnext;
    
    // helpful variables
    Chol<T> Qxx_chol;
    Chol<T> Quu_chol;
    DMat<T> Ixx;
    DMat<T> Iuu;
    DVec<T> _x0; // initial condition of current phase

public:
    size_t _xsize, _usize, _ysize;

};

#endif //SINGLEPHASE_ABSTRACT_H