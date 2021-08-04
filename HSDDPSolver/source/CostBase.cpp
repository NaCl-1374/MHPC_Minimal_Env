#include "CostBase.h"

template <typename TH, size_t XSIZE, size_t USIZE, size_t YSIZE>
void Cost<TH, XSIZE, USIZE, YSIZE>::running_cost(ModelState<TH, XSIZE, USIZE, YSIZE> &mstate,
                                                ModelState<TH, XSIZE, USIZE, YSIZE> &mstate_ref,
                                                RCostStruct<TH, XSIZE, USIZE, YSIZE> &rcost, int mode)
{
    assert(("Actual state size and reference state size not equal", mstate.x.size() == mstate_ref.x.size()));
    assert(("Actual control size and reference control size not equal", mstate.u.size() == mstate_ref.u.size()));
    assert(("Actual output size and reference output size not equal", mstate.y.size() == mstate_ref.y.size()));

    rcost.l = (mstate.x - mstate_ref.x).transpose() * this->_Q[mode - 1] * (mstate.x - mstate_ref.x);
    rcost.l += (mstate.u - mstate_ref.u).transpose() * this->_R[mode - 1] * (mstate.u - mstate_ref.u);
    rcost.l += (mstate.y - mstate_ref.y).transpose() * this->_S[mode - 1] * (mstate.y - mstate_ref.y);
    rcost.l *= _dt;
}

template <typename TH, size_t XSIZE, size_t USIZE, size_t YSIZE>
void Cost<TH, XSIZE, USIZE, YSIZE>::running_cost_par(ModelState<TH, XSIZE, USIZE, YSIZE> &mstate,
                                                    ModelState<TH, XSIZE, USIZE, YSIZE> &mstate_ref,
                                                    RCostStruct<TH, XSIZE, USIZE, YSIZE> &rcost, int mode)
{
    assert(("Actual state size and reference state size not equal", mstate.x.size() == mstate_ref.x.size()));
    assert(("Actual control size and reference control size not equal", mstate.u.size() == mstate_ref.u.size()));
    assert(("Actual output size and reference output size not equal", mstate.y.size() == mstate_ref.y.size()));

    rcost.lx = 2 * _dt * this->_Q[mode - 1] * (mstate.x - mstate_ref.x);
    rcost.lu = 2 * _dt * this->_R[mode - 1] * (mstate.u - mstate_ref.u);
    rcost.ly = 2 * _dt * this->_S[mode - 1] * (mstate.y - mstate_ref.y);
    rcost.lxx = 2 * _dt * this->_Q[mode - 1];
    rcost.luu = 2 * _dt * this->_R[mode - 1];
    rcost.lux.setZero();
    rcost.lyy = 2 * _dt * this->_S[mode - 1];    
}

template <typename TH, size_t XSIZE, size_t USIZE, size_t YSIZE>
void Cost<TH, XSIZE, USIZE, YSIZE>::terminal_cost(ModelState<TH, XSIZE, USIZE, YSIZE> & mstate,
                                                 ModelState<TH, XSIZE, USIZE, YSIZE> & mstate_ref,
                                                 TCostStruct<TH, XSIZE> & tcost, int mode)
{
    assert(("Actual state size and reference state size not equal", mstate.x.size() == mstate_ref.x.size()));
    assert(("Actual control size and reference control size not equal", mstate.u.size() == mstate_ref.u.size()));
    assert(("Actual output size and reference output size not equal", mstate.y.size() == mstate_ref.y.size()));

    tcost.Phi = (mstate.x - mstate_ref.x).transpose() * this->_Qf[mode - 1] * (mstate.x - mstate_ref.x);
    tcost.Phi *= 0.5;
}

template <typename TH, size_t XSIZE, size_t USIZE, size_t YSIZE>
void Cost<TH, XSIZE, USIZE, YSIZE>::terminal_cost_par(ModelState<TH, XSIZE, USIZE, YSIZE> & mstate,
                                                     ModelState<TH, XSIZE, USIZE, YSIZE> & mstate_ref,
                                                     TCostStruct<TH, XSIZE> & tcost, int  mode)
{
    assert(("Actual state size and reference state size not equal", mstate.x.size() == mstate_ref.x.size()));
    assert(("Actual control size and reference control size not equal", mstate.u.size() == mstate_ref.u.size()));
    assert(("Actual output size and reference output size not equal", mstate.y.size() == mstate_ref.y.size()));

    tcost.Phix = this->_Qf[mode - 1] * (mstate.x - mstate_ref.x);
    tcost.Phixx = this->_Qf[mode - 1];
}

template class CostAbstract<double>;
template class Cost<double, xsize_WB, usize_WB, ysize_WB>;
template class Cost<double, xsize_FB, usize_FB, ysize_FB>;
