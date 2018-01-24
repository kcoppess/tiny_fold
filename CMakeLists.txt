cmake_minimum_required(VERSION 2.8.12)                                                                                                                                                                                                         
project(test)

# forcing  -fPIC so we can make a shared object
set(CMAKE_BUILD_TYPE Release)
set_property(GLOBAL PROPERTY POSITION_INDEPENDENT_CODE ON)
# pybind11 setup 

# Include path for Python header files
include_directories(/usr/include/python2.7)

find_library(PYBIND11 pybind11)
set(PYBIND11_CPP_STANDARD -std=c++11)
# requires pybind11 to be in local directory need to find workaround!
add_subdirectory(pybind11/)

add_library(tinyfold MODULE rna.cc)
target_link_libraries(tinyfold pybind11::module)
