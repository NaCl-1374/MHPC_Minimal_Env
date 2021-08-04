#include "MHPCCost.h"

template <typename TH>
WBCost<TH>::WBCost(TH dt) : Cost<TH, xsize_WB, usize_WB, ysize_WB>(dt, 4)
{
    this->_Q = new MatMN<TH, xsize_WB, xsize_WB>[4];
    this->_R = new MatMN<TH, usize_WB, usize_WB>[4];
    this->_S = new MatMN<TH, ysize_WB, ysize_WB>[4];
    this->_Qf = new MatMN<TH, xsize_WB, xsize_WB>[4];
    set_weighting_matrices();
}

template <typename TH>
FBCost<TH>::FBCost(TH dt) : Cost<TH, xsize_FB, usize_FB, ysize_FB>(dt, 4)
{
    this->_Q = new MatMN<TH, xsize_FB, xsize_FB>[4];
    this->_R = new MatMN<TH, usize_FB, usize_FB>[4];
    this->_S = new MatMN<TH, ysize_FB, ysize_FB>[4];
    this->_Qf = new MatMN<TH, xsize_FB, xsize_FB>[4];
    set_weighting_matrices();
}

template <typename TH>
void WBCost<TH>::set_weighting_matrices()
{
    // Set weighting matrices for WB running and terminal cost
    VecM<TH, xsize_WB> q, qf[4];
    VecM<TH, usize_WB> r[4];
    VecM<TH, ysize_WB> s[4];
    q << 0, 10, 5, 4, 4, 4, 4, 2, 1, .01, 6, 6, 6, 6;

    qf[0] << 0, 20, 8, 3, 3, 3, 3, 3, 2, 0.01, 5, 5, 0.01, 0.01;
    qf[1] << 0, 20, 8, 3, 3, 3, 3, 3, 2, 0.01, 5, 5, 5, 5;
    qf[2] << 0, 20, 8, 3, 3, 3, 3, 3, 2, 0.01, 0.01, 0.01, 5, 5;
    qf[3] << 0, 20, 8, 3, 3, 3, 3, 3, 2, 0.01, 5, 5, 5, 5;

    r[0] << 5, 5, 1, 1;
    r[1].setOnes();
    r[2] << 1, 1, 5, 5;
    r[3].setOnes();

    std::fill(&s[0], &s[3], VecM<TH, ysize_WB>(0, 0, 0, 0));
    s[0].tail(2) << 0.3, 0.3;
    s[2].head(2) << 0.15, 0.15;

    for (size_t i = 0; i < this->_n_modes; i++)
    {
        this->_Q[i] = 0.01 * q.asDiagonal();
        this->_R[i] = 0.5 * r[i].asDiagonal();
        this->_S[i] = s[i].asDiagonal();
        this->_Qf[i] = 100 * qf[i].asDiagonal();
    }
}

template <typename TH>
void FBCost<TH>::set_weighting_matrices()
{
    VecM<TH, xsize_FB> q, qf;
    VecM<TH, usize_FB> r[4];
    q.setZero();
    qf.setZero();
    std::memset(&r[0], 0, 4 * sizeof(VecM<TH, usize_FB>));
    q << 0, 10, 5, 2, 1, 0.01;
    qf << 1, 20, 8, 3, 1, 0.01;
    r[0].tail(2) << 0.01, 0.01;
    r[2].head(2) << 0.01, 0.01;

    for (size_t i = 0; i < 4; i++)
    {
        this->_Q[i] = 0.01 * q.asDiagonal();
        this->_R[i] = r[i].asDiagonal();
        this->_Qf[i] = 100 * qf.asDiagonal();
        this->_S[i].setZero();
    }
}

template class WBCost<double>;
template class FBCost<double>;
