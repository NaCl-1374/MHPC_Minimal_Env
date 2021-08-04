#ifndef PLANARQUADRUPED
#define PLANARQUADRUPED

#include <vector>
#include <eigen3/Eigen/Dense>
#include "MHPC_CPPTypes.h"
#include "PlanarRobot.h"

using namespace std;

namespace cheetah2D
{
    constexpr size_t num_act_joint = 4;
    constexpr size_t num_config = 7;
    constexpr size_t num_leg = 2;
}

namespace linkID2D
{
    constexpr size_t body = 0;
    constexpr size_t F_hip = 1;
    constexpr size_t F_knee = 2;
    constexpr size_t H_hip = 3;
    constexpr size_t H_knee = 4;
    constexpr size_t F_foot = 5;
    constexpr size_t H_foot = 6;
}

namespace legID
{
    constexpr size_t FLEG = 0;
    constexpr size_t HLEG = 1;
}

template <typename T>
class PlanarQuadruped : public RobotBase<T>
{
public:
    PlanarQuadruped(T dt = 0.001) : RobotBase<T>(dt,
                                                 2 * cheetah2D::num_config,
                                                 cheetah2D::num_act_joint,
                                                 2 * cheetah2D::num_leg,
                                                 ModelType::FULL)
    {
        _Ac.setZero();
        _Bc.setZero();
        _xdot.setZero();
    }

private:
    // Intermediate continuous-time variables
    MatMN<T, xsize_WB, xsize_WB> _Ac;
    MatMN<T, xsize_WB, usize_WB> _Bc;
    VecM<T, xsize_WB> _xdot;

protected:
    // resolve the problem of hidden overloaded functions
    using RobotBase<T>::dynamics;
    using RobotBase<T>::dynamics_par;
    using RobotBase<T>::resetmap;
    using RobotBase<T>::resetmap_par;
    using RobotBase<T>::plan_foothold;

public:
    void dynamics(VecM<T, xsize_WB> &x, VecM<T, usize_WB> &u, VecM<T, xsize_WB> &x_next, VecM<T, ysize_WB> &y, int mode) override;
    void resetmap(VecM<T, xsize_WB> &x, VecM<T, xsize_WB> &x_next, VecM<T, ysize_WB> &y, int mode) override;
    void dynamics_par(VecM<T, xsize_WB> &x, VecM<T, usize_WB> &u, MatMN<T, xsize_WB, xsize_WB> &A, MatMN<T, xsize_WB, usize_WB> &B,
                      MatMN<T, ysize_WB, xsize_WB> &C, MatMN<T, ysize_WB, usize_WB> &D, int mode) override;
    void resetmap_par(VecM<T, xsize_WB> &x, MatMN<T, xsize_WB, xsize_WB> &Px, int mode) override;
    void plan_foothold(DVec<T> &x, T stance_time, int mode) override {}
    void FootJacobian(VecM<T, xsize_WB> &x, MatMN<T,2,qsize_WB> &J, MatMN<T,2,qsize_WB> &Jd, int foot);
    void linkJacobian(VecM<T, xsize_WB> &x, MatMN<T,2,qsize_WB> &J, MatMN<T,2,qsize_WB> &Jd, size_t linkidx, VecM<T,2> contactLoc){}
    void build_quadruped();
    MatMN<T,4,4> get_homoTransformation(VecM<T, qsize_WB> &q, size_t linkidx);
    VecM<T,2> get_contact_position(VecM<T, qsize_WB> &q, size_t linkidx, VecM<T,2> &contactLoc);
    VecM<T,2> get_leg_ext_vec(VecM<T, qsize_WB> &q, int legidx);

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    T bodyMass, hipLinkMass, kneeLinkMass; // currently not used
    T bodyLength, hipLinkLength, kneeLinkLength;
    T bodyWidth, hipLinkWidth, kneeLinkWidth;
    VecM<T, 2> bodyCoM, hipLinkCoM, kneeLinkCoM;
    VecM<T, 2> hipLoc, kneeLoc;
};

#endif
