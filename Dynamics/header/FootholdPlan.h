#ifndef FOOTHOLD_PLANNER
#define FOOTHOLD_PLANNER
#include "MHPC_CPPTypes.h"
#include "PlanarQuadruped.h"

template <typename T>
class FootholdPlanner
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW    
private:
    T _velcmd; // desired forward speed
    T _ground_height;
    VecM<T, qsize_WB> _config; // generalized joints of full planar quadruped
    PlanarQuadruped<T> *_robot = nullptr;

public:
    FootholdPlanner(PlanarQuadruped<T> *PQuad, T velcmd = 0, T ground = -0.404)
    {
        _velcmd = velcmd;
        _robot = PQuad;
        _ground_height = ground;
        _config.setZero();
    }

    template <typename Derived>
    VecM<T, 4> get_foothold_location(const Eigen::DenseBase<Derived> &pos, T stance_time, int mode)
    {
        assert(pos.size() == 3); // if the pos dimension is not 3, throw an error
        VecM<T, 2> pos_hip(0, 0);
        VecM<T, 4> foothold(0, 0, 0, 0);
        _config.head(3) = pos;
        VecM<T, 2> contact_location = VecM<T, 2>::Zero();

        switch (mode)
        {
        case 1:
            pos_hip = _robot->get_contact_position(_config, linkID2D::H_hip, contact_location);
            foothold[2] = pos_hip[0] + _velcmd * stance_time / 2;
            foothold[3] = _ground_height;
            break;
        case 3:
            pos_hip = _robot->get_contact_position(_config, linkID2D::F_hip, contact_location);
            foothold[0] = pos_hip[0] + _velcmd * stance_time / 2;
            foothold[1] = _ground_height;
            break;
        }

        return foothold;
    }
};

// template class FootholdPlanner<float>;
#endif