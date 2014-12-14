#Define the project we are working on
project(mscomplex-tri)

# Look for various software we need
cmake_minimum_required(VERSION 2.6)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} " -frounding-math -std=c++0x")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

add_subdirectory(dmsc)

option(BUILD_VIEWER "build the mscomplex viewer" OFF)

if(BUILD_VIEWER)
add_subdirectory(utls)
add_subdirectory(viewer)
endif(BUILD_VIEWER)



