### Dependencies
* [Eigen3](https://eigen.tuxfamily.org/)

### Build Instructions
```
mkdir build
cd build
cmake ..
make 
```

### Run Instructions
./mhpc_ctrl

You could set up the number of whole-body phases and number of rigid-body phases with MHPC_UserParameter, and configure
HSDDP option such as convergence threshold, maximum number of iterations using HSDDP_OPTION.