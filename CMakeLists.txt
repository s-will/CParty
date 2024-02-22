# Ian Wark July 6 2017

cmake_minimum_required(VERSION 3.1)
project(CParty)

# enable c++11
# set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \ -O3")


# create simfold library
add_subdirectory(simfold)

set(SOURCE
    W_final.h W_final.cpp 
    pseudo_loop.h pseudo_loop.cpp 
    V_final.h V_final.cpp 
    VM_final.h VM_final.cpp
    h_common.h h_common.cpp 
    hfold_validation.h hfold_validation.cpp 
    hfold.h HFold.cpp
)

# pass SIMFOLD_HOME to binary
set(SIMFOLD_HOME ${CMAKE_SOURCE_DIR}/simfold)
add_definitions("-DSIMFOLD_HOME=\"${SIMFOLD_HOME}\"")

# set path for executable to be created
# set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
# create executable
add_executable(CParty ${SOURCE})
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "CParty")

# enable C++11
target_compile_features(CParty PRIVATE cxx_nullptr)

# link to simfold
target_link_libraries(CParty LINK_PUBLIC simfold)
