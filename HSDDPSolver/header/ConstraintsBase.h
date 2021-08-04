#ifndef CONSTRAINTSBASE_H
#define CONSTRAINTSBASE_H

#include "MHPC_CPPTypes.h"
#include "MHPC_CompoundTypes.h"
#include "PlanarRobot.h"

#include <vector>

template <typename TH>
class Constraint
{
protected:
    RobotBase<TH> *_model;

    std::vector<AL_REB_PARAMETER<TH>> _params;      // al and reb params for every phase

    std::vector<size_t> _num_tconstr, _num_pconstr; // number of constraints in every phase

public:
    Constraint(){}

    Constraint(RobotBase<TH> *model) : _model(model){}

    virtual ~ Constraint() = default;

    void get_AL_REB_PARAMS(AL_REB_PARAMETER<TH> &param, int mode) { param = _params[mode - 1]; }

    size_t get_num_pconstraint(int modeidx) // get # pconstraint for mode modeidx
    {
        assert((modeidx <= _num_pconstr.size() && modeidx > 0));
        return _num_pconstr[modeidx - 1];
    }

    size_t get_num_tconstraint(int modeidx) // get # tconstraint for mode modeidx
    {
        assert((modeidx <= _num_tconstr.size() && modeidx > 0));
        return _num_tconstr[modeidx - 1];
    }

    virtual void terminal_constraint(ModelState<TH, xsize_WB, usize_WB, ysize_WB> &, TConstrStruct<TH, xsize_WB>*, int) {}

    virtual void terminal_constraint(ModelState<TH, xsize_FB, usize_FB, ysize_FB> &, TConstrStruct<TH, xsize_FB>*, int) {}

    virtual void path_constraint(ModelState<TH, xsize_WB, usize_WB, ysize_WB> &,  IneqConstrStruct<TH, xsize_WB, usize_WB, ysize_WB>*, int) {}

    virtual void path_constraint(ModelState<TH, xsize_FB, usize_FB, ysize_FB> &,  IneqConstrStruct<TH, xsize_FB, usize_FB, ysize_FB>*, int) {}

    virtual void initialize_AL_REB_PARAMS(){};
};


# endif // CONSTRAINTSBASE_H