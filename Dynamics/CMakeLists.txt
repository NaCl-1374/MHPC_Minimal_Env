file(GLOB_RECURSE sources "source/*.cpp")

add_library(dynamics SHARED ${sources})

target_include_directories(dynamics PUBLIC 
                                    "header/"
                                    "../../../common/include" # needed orientation tools 
                                    "../Common/header"                                    
                                    "../CasadiGen/header")

target_link_libraries(dynamics casadigen)
