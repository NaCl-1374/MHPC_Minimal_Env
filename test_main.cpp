#include "MultiPhaseDDP.h"
#include "MHPC_CPPTypes.h"
#include "MHPC_CompoundTypes.h"
#include "PlanarFloatingBase.h"
#include "PlanarQuadruped.h"
#include "Gait.h"
#include "ReferenceGen.h"
#include "MHPCConstraints.h"
#include "MHPCCost.h"
#include "MHPCLocomotion.h"

int main(int argn, char *argv[])
{
	
    // Choose HSDDP options
    HSDDP_OPTION<double> option;
    option.ReB_active = 1;
    option.AL_active = 1;
    option.max_AL_iter = 2;
    option.max_DDP_iter = 3;

    // Instantiate MHPCLocomotion
    USRCMD usrcmd;
    usrcmd.vel = 1.5;
    usrcmd.height = 0;
    usrcmd.roll=0;
    usrcmd.pitch = 0;
    usrcmd.yaw = 0;
    MHPCUserParameters mhpcparams;   
    mhpcparams.usrcmd = &usrcmd;
    Gait gait; 
    MHPCLocomotion<double> locomotion(&mhpcparams, &gait, option);    
    locomotion.initialization();
    locomotion.solve_mhpc();
    locomotion.print_debugInfo();
}