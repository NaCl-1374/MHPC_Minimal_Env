#ifndef REFERENCEGEN_H
#define REFERENCEGEN_H

#include "MHPC_CPPTypes.h"
#include "MHPC_CompoundTypes.h"
#include "PlanarRobot.h"
#include <iostream>

template <typename T>
class ReferenceGen
{
private:
    int _n_phases, _n_WBphases, _n_FBphases;
    T _dt_WB, _dt_FB;
    DVec<int> _mode_seq;
    DVec<int> _N_TIMESTEPS;
    vector<ModelState<T, xsize_WB, usize_WB, ysize_WB> *> _ref_WB;
    vector<ModelState<T, xsize_FB, usize_FB, ysize_FB> *> _ref_FB;
    T _vel_cmd, _height_cmd;

    VecM<T, 4> _qjoint_bias, _qd_joint_bias; // joint and joint vel const reference at running cost stage
    VecM<T, xsize_WB> _xref_WB_term[4];      // terminal state reference for WB planning
    T **_pos_ref = nullptr;
    T _GRF;

public:
    ReferenceGen() {}
    ReferenceGen(int n_WBphases, int n_FBphases,
                 vector<ModelState<T, xsize_WB, usize_WB, ysize_WB> *> &ref_WB,
                 vector<ModelState<T, xsize_FB, usize_FB, ysize_FB> *> &ref_FB,
                 T dt_WB, T dt_FB,
                 DVec<int> mode_seq,
                 DVec<int> N_TIMESTEPS,
                 T vel_cmd, T height_cmd);

    void Initialization(int n_WBphases, int n_FBphases,
                        vector<ModelState<T, xsize_WB, usize_WB, ysize_WB> *> &ref_WB,
                        vector<ModelState<T, xsize_FB, usize_FB, ysize_FB> *> &ref_FB,
                        T dt_WB, T dt_FB,
                        DVec<int> mode_seq,
                        DVec<int> N_TIMESTEPS,
                        T vel_cmd, T height_cmd);

    void update_phase_config(DVec<int> mode_seq, DVec<int> N_TIMESTEPS){_mode_seq = mode_seq; _N_TIMESTEPS = N_TIMESTEPS;}  

    void update_data(vector<ModelState<T, xsize_WB, usize_WB, ysize_WB> *> &ref_WB,
                     vector<ModelState<T, xsize_FB, usize_FB, ysize_FB> *> &ref_FB)  
                     {
                         _ref_WB = ref_WB;
                         _ref_FB = ref_FB;
                     }                    

    template <typename Derived>
    void generate_ref(const Eigen::DenseBase<Derived> &x0)
    {
        int modeidx;

        calc_forward_ref_pos(x0);

        for (size_t pidx(0); pidx < _n_phases; pidx++)
        {
            modeidx = _mode_seq[pidx];
            for (size_t k(0); k < _N_TIMESTEPS[pidx] - 1; k++)
            {
                if (pidx < _n_WBphases)
                {
                    _ref_WB[pidx][k].x << _pos_ref[pidx][k], _height_cmd, 0, _qjoint_bias,
                        _vel_cmd, 0, 0, _qd_joint_bias;
                    _ref_WB[pidx][k].y << 0, _GRF, 0, _GRF;
                }
                else
                {
                    _ref_FB[pidx - _n_WBphases][k].x << _pos_ref[pidx][k], _height_cmd, 0,
                        _vel_cmd, 0, 0;
                    _ref_FB[pidx - _n_WBphases][k].u << 0, _GRF, 0, _GRF;
                }
            }

            /* set terminal state reference */
            if (pidx < _n_WBphases)
            {
                _ref_WB[pidx][_N_TIMESTEPS[pidx] - 1].x = _xref_WB_term[modeidx - 1];
                _ref_WB[pidx][_N_TIMESTEPS[pidx] - 1].x[0] = _pos_ref[pidx][_N_TIMESTEPS[pidx] - 1];
            }
            else
            {
                _ref_FB[pidx - _n_WBphases][_N_TIMESTEPS[pidx] - 1].x.tail(5) << _height_cmd, 0, _vel_cmd, 0, 0;
                _ref_FB[pidx - _n_WBphases][_N_TIMESTEPS[pidx] - 1].x[0] = _pos_ref[pidx][_N_TIMESTEPS[pidx] - 1];
            }
        }
    }

private:
    template <typename Derived>
    void calc_forward_ref_pos(const Eigen::DenseBase<Derived> &x0)
    {
        T dt = 0.0;
        _pos_ref[0][0] = x0[0];
        for (size_t pidx(0); pidx < _n_phases; pidx++)
        {
            dt = (pidx < _n_WBphases) ? (_dt_WB) : (_dt_FB);
            if (pidx > 0)
                _pos_ref[pidx][0] = _pos_ref[pidx - 1][_N_TIMESTEPS[pidx - 1] - 1];
            for (size_t k = 1; k < _N_TIMESTEPS[pidx]; k++)
            {
                _pos_ref[pidx][k] = _pos_ref[pidx][k - 1] + _vel_cmd * dt;
            }
        }
    }

    void free_memory()
    {
        // delete dynamically allocated memory
        if (nullptr != _pos_ref)
        {
            for (size_t pidx = 0; pidx < _n_phases; pidx++)
            {
                if (nullptr != _pos_ref[pidx])
                {
                    delete[] _pos_ref[pidx];
                }
            }
            delete[] _pos_ref;
            _pos_ref = nullptr;
        }
    }

public:
    ~ReferenceGen()
    {
        free_memory();
    }
};

#endif