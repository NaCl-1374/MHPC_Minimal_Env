#ifndef COSTBASE_H
#define COSTBASE_H

#include "MHPC_CPPTypes.h"
#include "MHPC_CompoundTypes.h"
#include "PlanarRobot.h"
#include <vector>

template <typename TH>
class CostAbstract
{
public:
    virtual ~ CostAbstract() = default;

    virtual void running_cost(ModelState<TH, xsize_WB, usize_WB, ysize_WB> &,
                              ModelState<TH, xsize_WB, usize_WB, ysize_WB> &,
                              RCostStruct<TH, xsize_WB, usize_WB, ysize_WB> &, int) {}

    virtual void running_cost(ModelState<TH, xsize_FB, usize_FB, ysize_FB> &,
                              ModelState<TH, xsize_FB, usize_FB, ysize_FB> &,
                              RCostStruct<TH, xsize_FB, usize_FB, ysize_FB> &, int) {}

    virtual void running_cost_par(ModelState<TH, xsize_WB, usize_WB, ysize_WB> &,
                                  ModelState<TH, xsize_WB, usize_WB, ysize_WB> &,
                                  RCostStruct<TH, xsize_WB, usize_WB, ysize_WB> &, int) {}

    virtual void running_cost_par(ModelState<TH, xsize_FB, usize_FB, ysize_FB> &,
                                  ModelState<TH, xsize_FB, usize_FB, ysize_FB> &,
                                  RCostStruct<TH, xsize_FB, usize_FB, ysize_FB> &, int) {}

    virtual void terminal_cost(ModelState<TH, xsize_WB, usize_WB, ysize_WB> &,
                               ModelState<TH, xsize_WB, usize_WB, ysize_WB> &,
                               TCostStruct<TH, xsize_WB> &, int) {}

    virtual void terminal_cost(ModelState<TH, xsize_FB, usize_FB, ysize_FB> &,
                               ModelState<TH, xsize_FB, usize_FB, ysize_FB> &,
                               TCostStruct<TH, xsize_FB> &, int) {}

    virtual void terminal_cost_par(ModelState<TH, xsize_WB, usize_WB, ysize_WB> &,
                                   ModelState<TH, xsize_WB, usize_WB, ysize_WB> &,
                                   TCostStruct<TH, xsize_WB> &, int) {}

    virtual void terminal_cost_par(ModelState<TH, xsize_FB, usize_FB, ysize_FB> &,
                                   ModelState<TH, xsize_FB, usize_FB, ysize_FB> &,
                                   TCostStruct<TH, xsize_FB> &, int) {}
};

template <typename TH, size_t XSIZE, size_t USIZE, size_t YSIZE>
class Cost : public CostAbstract<TH>
{
public:
    TH _dt;
    MatMN<TH, XSIZE, XSIZE> *_Q = nullptr, *_Qf = nullptr;
    MatMN<TH, USIZE, USIZE> *_R = nullptr;
    MatMN<TH, YSIZE, YSIZE> *_S = nullptr;
    int _n_modes;
    static constexpr size_t _xsize = XSIZE, _usize = USIZE, _ysize = YSIZE;

public:
    Cost(TH dt=0.001, int n_modes=4):_dt(dt), _n_modes(n_modes) {}

    virtual void set_weighting_matrices(){}

protected:
    using CostAbstract<TH>::running_cost; // resolve problem of hidden overloaded functions
    using CostAbstract<TH>::running_cost_par;
    using CostAbstract<TH>::terminal_cost;
    using CostAbstract<TH>::terminal_cost_par;

public:
    virtual ~Cost() = default;
    
    void running_cost(ModelState<TH, XSIZE, USIZE, YSIZE> &,
                      ModelState<TH, XSIZE, USIZE, YSIZE> &,
                      RCostStruct<TH, XSIZE, USIZE, YSIZE> &, int) override;

    void running_cost_par(ModelState<TH, XSIZE, USIZE, YSIZE> &,
                          ModelState<TH, XSIZE, USIZE, YSIZE> &,
                          RCostStruct<TH, XSIZE, USIZE, YSIZE> &, int) override;

    void terminal_cost(ModelState<TH, XSIZE, USIZE, YSIZE> &,
                       ModelState<TH, XSIZE, USIZE, YSIZE> &,
                       TCostStruct<TH, XSIZE> &, int) override;

    void terminal_cost_par(ModelState<TH, XSIZE, USIZE, YSIZE> &,
                           ModelState<TH, XSIZE, USIZE, YSIZE> &,
                           TCostStruct<TH, XSIZE> &, int) override;

};

#endif // COSTBASE_H