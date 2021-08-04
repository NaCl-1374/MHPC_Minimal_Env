#ifndef MHPC_LOCOMOTION_H
#define MHPC_LOCOMOTION_H
#include "MultiPhaseDDP.h"
#include "MHPC_CPPTypes.h"
#include "MHPC_CompoundTypes.h"
#include "PlanarFloatingBase.h"
#include "PlanarQuadruped.h"
#include "Gait.h"
#include "ReferenceGen.h"
#include "MHPCConstraints.h"
#include "MHPCCost.h"

template<typename TH>
class MHPCLocomotion:public MultiPhaseDDP<TH>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

public:
    MHPCLocomotion(MHPCUserParameters *Params, Gait*gait_in, HSDDP_OPTION<TH>option_in);                
    ~MHPCLocomotion();
    void warmstart() override; // override default warm_start
    void initialization();
    void solve_mhpc(); // modify the solve function in base class MultiphaseDDP
    void run();
    void print_debugInfo();

private:
    void memory_alloc();
    void memory_reset(); // set all data to zero
    void memory_free();
    void build_problem();
    void update_problem();


private:  
    // trajectory data
    ModelState<TH,xsize_WB, usize_WB, ysize_WB> ** ms_act_WB = nullptr;
    ModelState<TH,xsize_WB, usize_WB, ysize_WB> ** ms_nom_WB = nullptr;
    ModelState<TH,xsize_WB, usize_WB, ysize_WB> ** ms_ref_WB = nullptr;
    DynDerivative<TH, xsize_WB, usize_WB, ysize_WB> **dynpar_WB = nullptr; // dynamics partial along nominal trajectory
    RCostStruct<TH, xsize_WB, usize_WB, ysize_WB> **rcost_WB = nullptr;
    CostToGoStruct<TH, xsize_WB, usize_WB> **CTG_WB = nullptr;
    TCostStruct<TH, xsize_WB> *tcost_WB = nullptr;

    ModelState<TH,xsize_FB, usize_FB, ysize_FB> ** ms_act_FB = nullptr;
    ModelState<TH,xsize_FB, usize_FB, ysize_FB> ** ms_nom_FB = nullptr;
    ModelState<TH,xsize_FB, usize_FB, ysize_FB> ** ms_ref_FB = nullptr;
    DynDerivative<TH, xsize_FB, usize_FB, ysize_FB> **dynpar_FB = nullptr; // dynamics partial along nominal trajectory
    RCostStruct<TH, xsize_FB, usize_FB, ysize_FB> **rcost_FB = nullptr;
    CostToGoStruct<TH, xsize_FB, usize_FB> **CTG_FB = nullptr;
    TCostStruct<TH, xsize_FB> *tcost_FB = nullptr;   

    ModelState<TH,xsize_WB, usize_WB, ysize_WB> *ms_exec = nullptr; // solution bag for execution horizon
    CostToGoStruct<TH, xsize_WB, usize_WB> *CTG_exec = nullptr;

    // configuration data
    int n_wbphase, n_fbphase;
    double dt_wb, dt_fb;
    int cmode;
    DVec<int> mode_seq;
    DVec<float> timings;        // time duration of each phase
    DVec<int> N_TIMESTEPS;      // num of time steps of each phase    
    std::vector<int> pidx_WB;   // phase index for whole-body
    std::vector<int> pidx_FB;   // phase index for floating base

    // default initial condition for planar quadruped
    VecM<TH, 7> q0, qd0;
    VecM<TH, 14> x0; 


private:
    MHPCUserParameters *userParams = nullptr;
    PlanarQuadruped<TH> *wbmodel = nullptr;
    PlanarFloatingBase<TH> *fbmodel = nullptr;
    ReferenceGen<TH> refGen;
    FBConstraint<TH> *fbconstrGen = nullptr;
    WBConstraint<TH> *wbconstrGen= nullptr;
    WBCost<TH> *wbcostGen = nullptr;
    FBCost<TH> *fbcostGen = nullptr;
    FootholdPlanner<TH> * footholdplanner = nullptr;
    Gait * gait = nullptr;
};




#endif //MHPC_LOCOMOTION_H