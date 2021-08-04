# ifndef MHPC_CONSTRAINTS_H
# define MHPC_CONSTRAINTS_H

#include "MHPC_CPPTypes.h"
#include "MHPC_CompoundTypes.h"
#include "ConstraintsBase.h"
#include "PlanarRobot.h"

template <typename TH>
class FBConstraint : public Constraint<TH>
{
protected:// resolve problem of hidden overloaded functions
    using Constraint<TH>::terminal_constraint;
    using Constraint<TH>::path_constraint;

public:
    FBConstraint(RobotBase<TH> *model);

    void terminal_constraint(ModelState<TH, xsize_FB, usize_FB, ysize_FB> &, TConstrStruct<TH, xsize_FB>*, int) override {}

    void path_constraint(ModelState<TH, xsize_FB, usize_FB, ysize_FB> &, IneqConstrStruct<TH, xsize_FB, usize_FB, ysize_FB>*, int) override {}
};

template <typename TH>
class WBConstraint : public Constraint<TH>
{
public:
    size_t _num_torque_limit;
    size_t _num_joint_limit;
    size_t _num_GRF_constraint;
    DMat<TH> C_torque, C_joint, C_grf;
    DVec<TH> b_torque, b_joint, b_grf;
    float  _friccoeff = 0.5; // static friction coefficient

protected: // resolve problem of hidden overloaded functions
    using Constraint<TH>::terminal_constraint;
    using Constraint<TH>::path_constraint;

public:
    WBConstraint(RobotBase<TH> *model);

    void terminal_constraint(ModelState<TH, xsize_WB, usize_WB, ysize_WB> &, TConstrStruct<TH, xsize_WB>*, int) override;

    void path_constraint(ModelState<TH, xsize_WB, usize_WB, ysize_WB> &, IneqConstrStruct<TH, xsize_WB, usize_WB, ysize_WB>*, int) override;

    void torque_limit(ModelState<TH, xsize_WB, usize_WB, ysize_WB> &, IneqConstrStruct<TH, xsize_WB, usize_WB, ysize_WB>*, int); // torque limit constraint

    void GRF_constraint(ModelState<TH, xsize_WB, usize_WB, ysize_WB> &, IneqConstrStruct<TH, xsize_WB, usize_WB, ysize_WB>*, int); // non-negative normal GRF and friction constraint

    void joint_limit(ModelState<TH, xsize_WB, usize_WB, ysize_WB> &, IneqConstrStruct<TH, xsize_WB, usize_WB, ysize_WB>*, int); // joint limit constraint

    void initialize_AL_REB_PARAMS();
};

# endif // MHPC_CONSTRAINTS_H