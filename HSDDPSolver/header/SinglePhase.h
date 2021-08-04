#ifndef SINGLEPHASE_HSDDP
#define SINGLEPHASE_HSDDP

#include <string>
#include "SinglePhaseAbstract.h"
using std::string;

template<typename> class MultiPhaseDDP; //forward declaration of MultiPhaseDDP class template

template<typename T, size_t XSIZE, size_t USIZE, size_t YSIZE>
class SinglePhase: public SinglePhaseAbstract<T>
{
private:
    friend class MultiPhaseDDP<T>;
    ModelState<T, XSIZE, USIZE, YSIZE> *_ms_act = nullptr;       // actual trajectory
    ModelState<T, XSIZE, USIZE, YSIZE> *_ms_nom = nullptr;       // norminal trajectory
    DynDerivative<T, XSIZE, USIZE, YSIZE> *_dynpar = nullptr; // dynamics partial along nominal trajectory
    RCostStruct<T, XSIZE, USIZE, YSIZE> *_rcost = nullptr;
    CostToGoStruct<T, XSIZE, USIZE> *_CTG = nullptr;
    TCostStruct<T, XSIZE> *_tcost = nullptr;        
    ModelState<T, XSIZE, USIZE, YSIZE> *_ms_ref = nullptr; // reference trajectory

    TConstrStruct<T, XSIZE>* _tconstr = nullptr;    
    IneqConstrStruct<T, XSIZE, USIZE, YSIZE>* _pconstr = nullptr; 
    DVec<T> B, Bz, Bzz;
    size_t _num_pconstr=0, _num_tconstr=0;
    AL_REB_PARAMETER<T> _param; // AL and ReB parameter for current phase

protected: // resolve problem of hidden overloadded function
    using SinglePhaseAbstract<T>::set_data;
    using SinglePhaseAbstract<T>::set_reference;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    SinglePhase():SinglePhaseAbstract<T>(XSIZE, USIZE, YSIZE){}
    
    SinglePhase(RobotBase<T> *, CostAbstract<T> *, Constraint<T> *, HSDDP_OPTION<T>*, T, size_t, size_t, size_t);

    ~SinglePhase();

    void initialization();

    void set_data(ModelState<T, XSIZE, USIZE, YSIZE>*, 
                  ModelState<T, XSIZE, USIZE, YSIZE>*,
                  DynDerivative<T, XSIZE, USIZE, YSIZE> *,
                  RCostStruct<T, XSIZE, USIZE, YSIZE> *,
                  CostToGoStruct<T, XSIZE, USIZE> *,
                  TCostStruct<T, XSIZE> *,
                  ModelState<T, XSIZE, USIZE, YSIZE> *) override;

    // void set_reference(ModelState<T,XSIZE,USIZE,YSIZE> *ms_ref) override { _ms_ref = ms_ref; }

    void forward_sweep(T eps) override;

    void forward_sweep_dynamics_only(T eps) override; // only compute dynamics propagation

    void forward_sweep_partials_only() override; // only compute dynamics partials

    bool backward_sweep(T regularization) override;

    void get_CTG_currentPhase(T &dV, DVec<T> &G, DMat<T> &H) override;

    DVec<T> get_terminal_state() override {return _ms_act[this->_N_TIMESTEPS-1].x;}

    virtual void* get_nominal_ms_ptr() override;

    virtual void* get_CTG_info_ptr() override;

    void update_nominal_trajectory() override ;

    void udpate_AL_ReB_Param() override;

    void empty_bag() override;

    void set_data_to_zero();

    void initialize_AL_ReB_Params() override;

    T get_tconstr_violation() 
    {
        T sum(0);
        for (size_t i = 0; i < _num_tconstr; i++)
        {
            sum += _tconstr[i].get_violation_normsquare();
        }
        
        return sum;
    }

private:
    void reduced_barrier(IneqConstrStruct<T, XSIZE, USIZE, YSIZE> *Ineq, DVec<T> &delta, DVec<T> &B, DVec<T> &Bz, DVec<T> &Bzz); // z: vector of inequality constraint
                                                                                            // B: vector of reduced barrier B
                                                                                            // Bz: first-order derivative
    void update_running_cost_with_pconstr(RCostStruct<T, XSIZE, USIZE, YSIZE> &Rcost, 
                                          IneqConstrStruct<T, XSIZE, USIZE, YSIZE>* Ineq, 
                                          DVec<T> &delta, DVec<T> &eps,  int calcflag);

    void update_running_cost_with_smooth();

    void update_terminal_cost_with_tconstr(TCostStruct<T,XSIZE> &Tcost, 
                                           TConstrStruct<T,XSIZE>* Tconstr, 
                                           T sigma, DVec<T> &lambda, int calcflag);
};

#endif // SINGLEPHASE_HSDDP