#include "ReferenceGen.h"
#include <iostream>

template <typename T>
ReferenceGen<T>::ReferenceGen(int n_WBphases, int n_FBphases,
                              vector<ModelState<T, xsize_WB, usize_WB, ysize_WB> *> &ref_WB,
                              vector<ModelState<T, xsize_FB, usize_FB, ysize_FB> *> &ref_FB,
                              T dt_WB, T dt_FB,
                              DVec<int> mode_seq,
                              DVec<int> N_TIMESTEPS,
                              T vel_cmd, T height_cmd)
{
    Initialization(n_WBphases, n_FBphases,
                   ref_WB,
                   ref_FB,
                   dt_WB, dt_FB,
                   mode_seq,
                   N_TIMESTEPS,
                   vel_cmd, height_cmd);
}

template <typename T>
void ReferenceGen<T>::Initialization(int n_WBphases, int n_FBphases,
                                     vector<ModelState<T, xsize_WB, usize_WB, ysize_WB> *> &ref_WB,
                                     vector<ModelState<T, xsize_FB, usize_FB, ysize_FB> *> &ref_FB,
                                     T dt_WB, T dt_FB,
                                     DVec<int> mode_seq,
                                     DVec<int> N_TIMESTEPS,
                                     T vel_cmd, T height_cmd)
{
    _n_WBphases = n_WBphases;
    _n_FBphases = n_FBphases;
    _ref_WB = ref_WB;
    _ref_FB = ref_FB;
    _dt_WB = dt_WB;
    _dt_FB = dt_FB;
    _n_phases = n_WBphases + n_FBphases;
    _mode_seq = mode_seq;
    _N_TIMESTEPS = N_TIMESTEPS;
    _vel_cmd = vel_cmd;
    _height_cmd = height_cmd;
    _GRF = 8.252 * 9.81;
    std::memset(&_xref_WB_term[0], 0, 4 * sizeof(VecM<T, xsize_WB>));

    _xref_WB_term[0] << 0, -0.1432, -PI / 25, 0.35 * PI, -0.65 * PI, 0.35 * PI, -0.6 * PI,
        _vel_cmd, 1, 0, 0, 0, 0, 0;
    _xref_WB_term[1] << 0, -0.1418, PI / 35, 0.2 * PI, -0.58 * PI, 0.25 * PI, -0.7 * PI,
        _vel_cmd, -1, 0, 0, 0, 0, 0;
    _xref_WB_term[2] << 0, -0.1325, -PI / 40, 0.33 * PI, -0.48 * PI, 0.33 * PI, -0.75 * PI,
        _vel_cmd, 1, 0, 0, 0, 0, 0;
    _xref_WB_term[3] << 0, -0.1490, -PI / 25, 0.35 * PI, -0.7 * PI, 0.25 * PI, -0.60 * PI,
        _vel_cmd, -1, 0, 0, 0, 0, 0;

    _qjoint_bias << 0.3 * PI, -0.7 * PI, 0.3 * PI, -0.7 * PI;
    _qd_joint_bias.setZero();

    // allocate memory for forward position reference trajectory
    if (nullptr != _pos_ref)
        free_memory();

    _pos_ref = new T *[_n_phases];
    for (size_t pidx(0); pidx < _n_phases; pidx++)
    {
        _pos_ref[pidx] = new T[N_TIMESTEPS_MAX];
    }
}


template class ReferenceGen<double>;