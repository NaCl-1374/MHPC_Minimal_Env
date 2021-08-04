#include "MHPCLocomotion.h"
#include "boundingPDControl.h"
#include <unistd.h>
#include <iostream>
#include <fstream>

template<typename TH>
MHPCLocomotion<TH>::MHPCLocomotion(MHPCUserParameters*params, Gait * gait_in, HSDDP_OPTION<TH>option_in):                
                MultiPhaseDDP<TH>(params->n_wbphase+params->n_fbphase, 
                                 option_in),
                                 dt_wb(params->dt_wb),
                                 dt_fb(params->dt_fb),
                                 userParams(params),
                                 n_wbphase(params->n_wbphase),
                                 n_fbphase(params->n_fbphase),
                                 cmode(params->cmode)  ,
                                 gait(gait_in)                                   
{
    mode_seq.setOnes(this->_n_phases);
    timings.setZero(this->_n_phases);
    N_TIMESTEPS.setZero(this->_n_phases);

    wbmodel = new PlanarQuadruped<TH>(dt_wb);
    wbmodel->build_quadruped();
    footholdplanner = new FootholdPlanner<TH>(wbmodel, 1.5, -0.404);
    fbmodel = new PlanarFloatingBase<TH>(footholdplanner, dt_fb);
    wbcostGen = new WBCost<TH>(dt_wb);
    fbcostGen = new FBCost<TH>(dt_fb);
    wbconstrGen = new WBConstraint<TH>(wbmodel);
    fbconstrGen = new FBConstraint<TH>(fbmodel);

    this->_stateProj.setZero(xsize_FB, xsize_WB);
    this->_stateProj.block(0, 0, 3, 3) = DMat<TH>::Identity(3, 3);
    this->_stateProj.block(3, 7, 3, 3) = DMat<TH>::Identity(3, 3);

    /* set default initial condition (back stance phase)*/
    q0 << 0.0927, -0.1093, -0.1542, 1.0957, -2.2033, 0.9742, -1.7098;
    qd0 << 0.9011, 0.2756, 0.7333, 0.0446, 0.0009, 1.3219, 2.7346;
    x0 << q0, qd0;

    /* initialize class members */
    memory_alloc(); // allocate memory for datafile(GLOB_RECURSE test_sources "test/test_*.cpp")             # test cpp files
}


template<typename TH>
void MHPCLocomotion<TH>::initialization()
{
    this->set_initial_condition(x0);
    memory_reset();
    build_problem();
    warmstart();
}


/*
  @brief: build multi-phase problem
  @params:
          cmode: current mode index
          usrcmd: user command
*/
template<typename TH>
void MHPCLocomotion<TH>::build_problem()
{
    mode_seq = gait->get_mode_seq(cmode, this->_n_phases); // compute mode index for each phase
    timings = gait->get_timings(mode_seq);                 // compute timings for each phase

    vector<ModelState<TH, xsize_WB, usize_WB, ysize_WB> *> ref_WB_vec;
    vector<ModelState<TH, xsize_FB, usize_FB, ysize_FB> *> ref_FB_vec;

    /* configure single-phase problem */
    for (size_t pidx = 0; pidx < this->_n_phases; pidx++)
    {
        if (pidx < n_wbphase) // If full-model whole body phase
        {
            N_TIMESTEPS[pidx] = round(timings[pidx] / dt_wb);
            this->_phases[pidx]->set_models(wbmodel, wbcostGen, wbconstrGen, &this->_option, dt_wb);
            this->_phases[pidx]->set_phase_config(mode_seq[pidx], pidx, N_TIMESTEPS[pidx]);
            this->_phases[pidx]->set_data(ms_act_WB[pidx], ms_nom_WB[pidx], dynpar_WB[pidx], rcost_WB[pidx],
                                          CTG_WB[pidx], &tcost_WB[pidx], ms_ref_WB[pidx]);
            ref_WB_vec.push_back(ms_ref_WB[pidx]);
        }
        else // If simple-model rigid body phase
        {
            N_TIMESTEPS[pidx] = round(timings[pidx] / dt_fb);
            this->_phases[pidx]->set_models(fbmodel, fbcostGen, fbconstrGen, &this->_option, dt_fb);
            this->_phases[pidx]->set_phase_config(mode_seq[pidx], pidx, N_TIMESTEPS[pidx]);
            this->_phases[pidx]->set_data(ms_act_FB[pidx - n_wbphase], ms_nom_FB[pidx - n_wbphase], dynpar_FB[pidx - n_wbphase], rcost_FB[pidx - n_wbphase],
                                          CTG_FB[pidx - n_wbphase], &tcost_FB[pidx - n_wbphase], ms_ref_FB[pidx - n_wbphase]);
            ref_FB_vec.push_back(ms_ref_FB[pidx - n_wbphase]);
        }
        (pidx < this->_n_phases - 1) ? this->_phases[pidx]->set_next_phase(this->_phases[pidx + 1]) // assign next phase
                                     : this->_phases[pidx]->set_next_phase(nullptr);
        this->_phases[pidx]->initialization();
    }

    /* Initialize reference generator */
    refGen.Initialization(n_wbphase, n_fbphase,
                           ref_WB_vec, ref_FB_vec,
                           dt_wb, dt_fb,
                           mode_seq, N_TIMESTEPS, userParams->usrcmd->vel, userParams->usrcmd->height);
    /* computes tracking reference for multi-phase problem */
    refGen.generate_ref(this->_x0);
}

template<typename TH>
void MHPCLocomotion<TH>::update_problem()
{
    vector<ModelState<TH, xsize_WB, usize_WB, ysize_WB> *> ref_WB_vec;
    vector<ModelState<TH, xsize_FB, usize_FB, ysize_FB> *> ref_FB_vec;

    int pidx_temp = pidx_WB.front();
    pidx_WB.erase(pidx_WB.begin());
    pidx_WB.push_back(pidx_temp);
    if (n_fbphase > 0)
    {
        pidx_temp = pidx_FB.front();
        pidx_FB.erase(pidx_FB.begin());
        pidx_FB.push_back(pidx_temp);
    }

    cmode = gait->get_next_mode(cmode);
    mode_seq = gait->get_mode_seq(cmode, this->_n_phases); // update mode sequence
    timings = gait->get_timings(mode_seq);                 // update phase timing

    /* update mode index tied to each single phase */
    for (size_t idx = 0; idx < this->_n_phases; idx++)
    {
        N_TIMESTEPS[idx] = round(timings[idx] / this->_phases[idx]->_dt);
        this->_phases[idx]->set_phase_config(mode_seq[idx], idx, N_TIMESTEPS[idx]);
    }

    /* update data tied to each single phase */
    int pidx;
    for (size_t idx = 0; idx < this->_n_phases; idx++)
    {
        if (idx < n_wbphase) // If full-model whole body phase
        {
            pidx = pidx_WB[idx];
            this->_phases[idx]->set_data(ms_act_WB[pidx], ms_nom_WB[pidx], dynpar_WB[pidx], rcost_WB[pidx],
                                         CTG_WB[pidx], &tcost_WB[pidx], ms_ref_WB[pidx]);
            ref_WB_vec.push_back(ms_ref_WB[pidx]);                                         
        }
        else // If simple-model rigid body phase
        {
            pidx = pidx_FB[idx -n_wbphase];
            this->_phases[idx]->set_data(ms_act_FB[pidx], ms_nom_FB[pidx], dynpar_FB[pidx], rcost_FB[pidx],
                                          CTG_FB[pidx], &tcost_FB[pidx], ms_ref_FB[pidx]);
            ref_FB_vec.push_back(ms_ref_FB[pidx]);                                
        }
        this->_phases[idx]->initialization();
    }

    /* Initialize reference generator */
    refGen.update_phase_config(mode_seq, N_TIMESTEPS);
    refGen.update_data(ref_WB_vec, ref_FB_vec);
    refGen.generate_ref(this->_x0);
}



/*
    Reimplements the solve() function as defined in base class MultiPhaseDDP by copying the solution 
    to the execution horizon 
*/
template<typename TH>
void MHPCLocomotion<TH>::solve_mhpc()
{
#ifdef DEBUG
    printf("\x1b[32m---------------------------------\n");
    printf("          Solving MHPC  \n");
    printf("---------------------------------\x1b[0m \n");
#endif
    this->solve();
    assert(("Whole-body phase needs to be at least one", n_wbphase>=1));

    /* copy nominal state and feedback information to the execution horizon from the first phase of the solution bag*/
    std::memcpy(ms_exec, 
               (ModelState<TH, xsize_WB, usize_WB, ysize_WB>*)this->_phases[0]->get_nominal_ms_ptr(), 
               sizeof(ModelState<TH, xsize_WB, usize_WB, ysize_WB>)*this->_phases[0]->_N_TIMESTEPS);
    std::memcpy(CTG_exec, 
               (CostToGoStruct<TH, xsize_WB, usize_WB>*)this->_phases[0]->get_CTG_info_ptr(), 
               sizeof(CostToGoStruct<TH, xsize_WB, usize_WB>)*this->_phases[0]->_N_TIMESTEPS);
            
    /* copy the second phase if any*/
    if (n_wbphase > 1)
    {
        std::memcpy(ms_exec + this->_phases[0]->_N_TIMESTEPS, 
                   (ModelState<TH, xsize_WB, usize_WB, ysize_WB>*)this->_phases[1]->get_nominal_ms_ptr(), 
                   sizeof(ModelState<TH, xsize_WB, usize_WB, ysize_WB>)*this->_phases[1]->_N_TIMESTEPS);
        std::memcpy(CTG_exec + this->_phases[0]->_N_TIMESTEPS, 
                   (CostToGoStruct<TH, xsize_WB, usize_WB>*)this->_phases[1]->get_CTG_info_ptr(), 
                   sizeof(CostToGoStruct<TH, xsize_WB, usize_WB>)*this->_phases[1]->_N_TIMESTEPS);
    } 
}



template<typename TH>
void MHPCLocomotion<TH>::warmstart()
{
    this->empty_bag();
    DVec<TH> x0_phase = this->_x0;
    for (int idx = 0; idx < n_wbphase; idx++)
    {
        int pidx = pidx_WB[idx];
        ms_act_WB[pidx][0].x = x0_phase;
        bounding_PDcontrol(wbmodel, ms_act_WB[pidx], mode_seq[idx], N_TIMESTEPS[idx]);
        if (idx + 1 < n_wbphase)
        {
            x0_phase = this->phase_transition(this->_phases[idx], this->_phases[idx + 1]);
        }
    }
    this->update_nominal_trajectory();
}

template<typename TH>
void MHPCLocomotion<TH>::memory_alloc()
{
    ms_act_WB = new ModelState<TH, xsize_WB, usize_WB, ysize_WB> *[n_wbphase];
    ms_nom_WB = new ModelState<TH, xsize_WB, usize_WB, ysize_WB> *[n_wbphase];
    ms_ref_WB = new ModelState<TH, xsize_WB, usize_WB, ysize_WB> *[n_wbphase];
    dynpar_WB = new DynDerivative<TH, xsize_WB, usize_WB, ysize_WB> *[n_wbphase];
    CTG_WB = new CostToGoStruct<TH, xsize_WB, usize_WB> *[n_wbphase];
    rcost_WB = new RCostStruct<TH, xsize_WB, usize_WB, ysize_WB> *[n_wbphase];
    tcost_WB = new TCostStruct<TH, xsize_WB>[n_wbphase];

    ms_act_FB = new ModelState<TH, xsize_FB, usize_FB, ysize_FB> *[n_fbphase];
    ms_nom_FB = new ModelState<TH, xsize_FB, usize_FB, ysize_FB> *[n_fbphase];
    ms_ref_FB = new ModelState<TH, xsize_FB, usize_FB, ysize_FB> *[n_fbphase];
    dynpar_FB = new DynDerivative<TH, xsize_FB, usize_FB, ysize_FB> *[n_fbphase];
    CTG_FB = new CostToGoStruct<TH, xsize_FB, usize_FB> *[n_fbphase];
    rcost_FB = new RCostStruct<TH, xsize_FB, usize_FB, ysize_FB> *[n_fbphase];
    tcost_FB = new TCostStruct<TH, xsize_FB>[n_fbphase];

    for (size_t idx = 0; idx < n_wbphase; idx++)
    {
        pidx_WB.push_back(idx);
        this->_phases[idx] = new SinglePhase<TH, xsize_WB, usize_WB, ysize_WB>;
        ms_act_WB[idx] = new ModelState<TH, xsize_WB, usize_WB, ysize_WB>[N_TIMESTEPS_MAX];
        ms_nom_WB[idx] = new ModelState<TH, xsize_WB, usize_WB, ysize_WB>[N_TIMESTEPS_MAX];
        ms_ref_WB[idx] = new ModelState<TH, xsize_WB, usize_WB, ysize_WB>[N_TIMESTEPS_MAX];
        dynpar_WB[idx] = new DynDerivative<TH, xsize_WB, usize_WB, ysize_WB>[N_TIMESTEPS_MAX];
        CTG_WB[idx] = new CostToGoStruct<TH, xsize_WB, usize_WB>[N_TIMESTEPS_MAX];
        rcost_WB[idx] = new RCostStruct<TH, xsize_WB, usize_WB, ysize_WB>[N_TIMESTEPS_MAX];
    }

    for (size_t idx = 0; idx < n_fbphase; idx++)
    {
        pidx_FB.push_back(idx);
        this->_phases[idx + n_wbphase] = new SinglePhase<TH, xsize_FB, usize_FB, ysize_FB>;
        ms_act_FB[idx] = new ModelState<TH, xsize_FB, usize_FB, ysize_FB>[N_TIMESTEPS_MAX];
        ms_nom_FB[idx] = new ModelState<TH, xsize_FB, usize_FB, ysize_FB>[N_TIMESTEPS_MAX];
        ms_ref_FB[idx] = new ModelState<TH, xsize_FB, usize_FB, ysize_FB>[N_TIMESTEPS_MAX];
        dynpar_FB[idx] = new DynDerivative<TH, xsize_FB, usize_FB, ysize_FB>[N_TIMESTEPS_MAX];
        CTG_FB[idx] = new CostToGoStruct<TH, xsize_FB, usize_FB>[N_TIMESTEPS_MAX];
        rcost_FB[idx] = new RCostStruct<TH, xsize_FB, usize_FB, ysize_FB>[N_TIMESTEPS_MAX];
    }
    ms_exec = new ModelState<TH, xsize_WB, usize_WB, ysize_WB>[2*N_TIMESTEPS_MAX];
    CTG_exec = new CostToGoStruct<TH, xsize_WB, usize_WB>[2*N_TIMESTEPS_MAX];
}


template<typename TH>
void MHPCLocomotion<TH>::memory_reset()
{
    for (size_t idx = 0; idx < n_wbphase; idx++)
    {
        std::memset(ms_act_WB[idx], 0, sizeof(ModelState<TH, xsize_WB, usize_WB, ysize_WB>) * N_TIMESTEPS_MAX);
        std::memset(ms_nom_WB[idx], 0, sizeof(ModelState<TH, xsize_WB, usize_WB, ysize_WB>) * N_TIMESTEPS_MAX);
        std::memset(ms_ref_WB[idx], 0, sizeof(ModelState<TH, xsize_WB, usize_WB, ysize_WB>) * N_TIMESTEPS_MAX);
        std::memset(dynpar_WB[idx], 0, sizeof(DynDerivative<TH, xsize_WB, usize_WB, ysize_WB>) * N_TIMESTEPS_MAX);
        std::memset(rcost_WB[idx], 0, sizeof(RCostStruct<TH, xsize_WB, usize_WB, ysize_WB>) * N_TIMESTEPS_MAX);
        std::memset(CTG_WB[idx], 0, sizeof(CostToGoStruct<TH, xsize_WB, usize_WB>) * N_TIMESTEPS_MAX);
        tcost_WB[idx].Zeros();
    }

    for (size_t idx = 0; idx < n_fbphase; idx++)
    {
        std::memset(ms_act_FB[idx], 0, sizeof(ModelState<TH, xsize_FB, usize_FB, ysize_FB>) * N_TIMESTEPS_MAX);
        std::memset(ms_nom_FB[idx], 0, sizeof(ModelState<TH, xsize_FB, usize_FB, ysize_FB>) * N_TIMESTEPS_MAX);
        std::memset(ms_ref_FB[idx], 0, sizeof(ModelState<TH, xsize_FB, usize_FB, ysize_FB>) * N_TIMESTEPS_MAX);
        std::memset(dynpar_FB[idx], 0, sizeof(DynDerivative<TH, xsize_FB, usize_FB, ysize_FB>) * N_TIMESTEPS_MAX);
        std::memset(rcost_FB[idx], 0, sizeof(RCostStruct<TH, xsize_FB, usize_FB, ysize_FB>) * N_TIMESTEPS_MAX);
        std::memset(CTG_FB[idx], 0, sizeof(CostToGoStruct<TH, xsize_FB, usize_FB>) * N_TIMESTEPS_MAX);
        tcost_FB[idx].Zeros();
    }
    std::memset(ms_exec, 0, sizeof(ModelState<TH, xsize_WB, usize_WB, ysize_WB>) * 2*N_TIMESTEPS_MAX);
    std::memset(CTG_exec, 0, sizeof(CostToGoStruct<TH, xsize_WB, usize_WB>) * 2*N_TIMESTEPS_MAX);    
}

template<typename TH>
void MHPCLocomotion<TH>::print_debugInfo()
{
    std::ofstream gradient_output("gradient.txt");
    std::ofstream state_output("state.txt");
    std::ofstream contrl_output("control.txt");
    std::ofstream cost_output("cost.txt");
    if (state_output.is_open())
    {
        printf("********** Write to file state.txt ************\n");
        for (int i = 0; i < n_wbphase; i++)
        {
            for (size_t k = 0; k < N_TIMESTEPS[i]; k++)
            {
                state_output << ms_nom_WB[i][k].x.transpose() << endl;
            }
        }
        for (int i = 0; i < n_fbphase; i++)
        {
            for (size_t k = 0; k < N_TIMESTEPS[i+n_wbphase]; k++)
            {
                state_output << ms_nom_FB[i][k].x.transpose() << endl;
            }
        }
    }
    state_output.close();
    if (contrl_output.is_open())
    {
        printf("********** Write to file control.txt ************\n");
        for (int i = 0; i < n_wbphase; i++)
        {
            for (size_t k = 0; k < N_TIMESTEPS[i]; k++)
            {
                contrl_output << ms_nom_WB[i][k].u.transpose() << endl;
            }
        }
        for (int i = 0; i < n_fbphase; i++)
        {
            for (size_t k = 0; k < N_TIMESTEPS[i+2]; k++)
            {
                contrl_output << ms_nom_FB[i][k].u.transpose() << endl;
            }
        }

    }
    contrl_output.close();
    if (gradient_output.is_open())
    {
        printf("********** Write to file gradient.txt ************\n");
        for (int i = 0; i < n_wbphase; i++)
        {
            for (size_t k = 0; k < N_TIMESTEPS[i]; k++)
            {
                gradient_output << CTG_WB[i][k].G.transpose() << endl;
            }
        }
        for (int i = 0; i < n_fbphase; i++)
        {
            for (size_t k = 0; k < N_TIMESTEPS[i+2]; k++)
            {
                gradient_output << CTG_FB[i][k].G.transpose() << endl;
            }
        }

    }
    gradient_output.close();
    if (cost_output.is_open())
    {
        printf("********** Write to file cost.txt ************\n");
        for (int i = 0; i < n_wbphase; i++)
        {
            for (size_t k = 0; k < N_TIMESTEPS[i]-1; k++)
            {
                cost_output << rcost_WB[i][k].lx.transpose() << endl;
            }
            cost_output << tcost_WB[i].Phix.transpose() << endl;
        }
        for (int i = 0; i < n_fbphase; i++)
        {
            for (size_t k = 0; k < N_TIMESTEPS[i+2]-1; k++)
            {
                cost_output << rcost_FB[i][k].lx.transpose() << endl;
            }
            cost_output << tcost_FB[i].Phix.transpose() << endl;
        }

    }
    cost_output.close();
}


template<typename TH>
MHPCLocomotion<TH>::~MHPCLocomotion()
{
    memory_free();
}

template<typename TH>
void MHPCLocomotion<TH>::memory_free()
{
    for (size_t i = 0; i < n_wbphase; i++)
    {
        delete[] ms_act_WB[i];
        delete[] ms_nom_WB[i];
        delete[] ms_ref_WB[i];
        delete[] dynpar_WB[i];
        delete[] CTG_WB[i];
        delete[] rcost_WB[i];
    }

    delete[] ms_act_WB;
    delete[] ms_nom_WB;
    delete[] ms_ref_WB;
    delete[] dynpar_WB;
    delete[] CTG_WB;
    delete[] rcost_WB;
    delete[] tcost_WB;

    ms_act_WB = nullptr;
    ms_nom_WB = nullptr;
    ms_ref_WB = nullptr;
    dynpar_WB = nullptr;
    CTG_WB = nullptr;
    rcost_WB = nullptr;
    tcost_WB = nullptr;

    delete [] ms_exec;
    delete [] CTG_exec;
    ms_exec = nullptr;
    CTG_exec = nullptr;

    for (size_t i = 0; i < n_fbphase; i++)
    {
        delete[] ms_act_FB[i];
        delete[] ms_nom_FB[i];
        delete[] ms_ref_FB[i];
        delete[] dynpar_FB[i];
        delete[] CTG_FB[i];
        delete[] rcost_FB[i];
    }

    delete[] ms_act_FB;
    delete[] ms_nom_FB;
    delete[] ms_ref_FB;
    delete[] dynpar_FB;
    delete[] CTG_FB;
    delete[] rcost_FB;
    delete[] tcost_FB;

    ms_act_FB = nullptr;
    ms_nom_FB = nullptr;
    ms_ref_FB = nullptr;
    dynpar_FB = nullptr;
    CTG_FB = nullptr;
    rcost_FB = nullptr;
    tcost_FB = nullptr;

    delete wbmodel;
    delete fbmodel;
    delete fbconstrGen;
    delete wbconstrGen;
    delete wbcostGen;
    delete fbcostGen;

    wbmodel = nullptr;
    fbmodel = nullptr;
    fbconstrGen = nullptr;
    wbconstrGen= nullptr;
    wbcostGen = nullptr;
    fbcostGen = nullptr;

    for (size_t i = 0; i < this->_n_phases; i++)
    {
        delete this->_phases[i];
        this->_phases[i] = nullptr;
    }


}

template class MHPCLocomotion<double>;
