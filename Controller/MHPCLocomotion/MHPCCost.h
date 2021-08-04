# ifndef MHPCCOST_H
# define MHPCCOST_H

#include "MHPC_CPPTypes.h"
#include "MHPC_CompoundTypes.h"
#include "CostBase.h"

template <typename TH>
class WBCost : public Cost<TH, xsize_WB, usize_WB, ysize_WB>
{
public:
    WBCost(TH dt);
    void set_weighting_matrices() override;
};

template <typename TH>
class FBCost : public Cost<TH, xsize_FB, usize_FB, ysize_FB>
{
public:
    FBCost(TH dt);
    void set_weighting_matrices() override;
};

#endif //MHPCCOST_H
