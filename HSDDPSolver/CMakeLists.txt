cmake_minimum_required(VERSION 3.10)

file(GLOB_RECURSE sources "source/*.cpp")

add_library(hsddp SHARED ${sources})

target_include_directories(hsddp PUBLIC
								 "../Common/header"
                                 "../Dynamics/header"
                                 "header")

