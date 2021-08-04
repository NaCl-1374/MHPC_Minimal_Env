cmake_minimum_required(VERSION 3.10)

project(MHPCCTRL)

# set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(DEBUG "Do not debug" OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}\
					 -O3 -no-pie -march=native -ggdb\
					 -Wall -Wextra -Werror\
					 -Wno-unused -Wno-error=unused\
					 -Wno-return-type -Wno-error=return-type\
					 -Wno-sign-compare -Wno-error=sign-compare\
					 -Wno-switch -Wno-error=switch\
					 -Wno-reorder -Wno-error=reorder")					 
					 
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}\
					 -O3 -march=native -std=gnu99 -ggdb\
					 -Wno-unused -Wno-error=unused\
					 -Wno-return-type -Wno-error=return-type\
					 -Wno-sign-compare -Wno-error=sign-compare\
					 -Wno-switch -Wno-error=switch")

include_directories("Common/header")
include_directories("Controller/Heuristics")
include_directories("Controller/MHPCLocomotion")

add_subdirectory(CasadiGen)
add_subdirectory(Dynamics)
add_subdirectory(HSDDPSolver)


file(GLOB_RECURSE sources "Controller/Heuristics/*.cpp"
						  "Controller/MHPCLocomotion/*.cpp")

add_executable(mhpc_ctrl test_main.cpp ${sources})
target_link_libraries(mhpc_ctrl casadigen dynamics hsddp)



