#ifndef BOUNDING_PD_CONTROLLER_H
#define BOUNDING_PD_CONTROLLER_H
#include "PlanarQuadruped.h"
#include "MHPC_CPPTypes.h"
#include "MHPC_CompoundTypes.h"
#include <vector>

template<typename T>
void bounding_PDcontrol(PlanarQuadruped<T>*, ModelState<T,xsize_WB,usize_WB,ysize_WB> *ms, int modeidx, size_t);

#endif