project(mscomplex-tri-core)

find_package(Boost COMPONENTS program_options serialization system REQUIRED)
find_package(Eigen3 REQUIRED)

cmake_minimum_required(VERSION 2.8.8)

include_directories(
  ${EIGEN3_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${Boost_INCLUDE_DIRS})
  

set(MSCOMPLEX_TRI_CORE_SRCS

    trimesh.h

    trimesh_dataset.h
    trimesh_dataset_ensure.h
    trimesh_dataset.cpp

    trimesh_mscomplex.h
    trimesh_mscomplex_ensure.h
    trimesh_mscomplex_simp.h
    trimesh_mscomplex.cpp
#    trimesh_mscomplex_vtk.cpp

    utl.h
    utl.cpp

    tri_edge.h
    tri_edge.cpp
    )

set_source_files_properties(${MSCOMPLEX_TRI_CORE_SRCS} PROPERTIES COMPILE_FLAGS -fpic)

add_library(mscomplex-tri-core OBJECT ${MSCOMPLEX_TRI_CORE_SRCS})

option(BUILD_CL_TOOL "command line tool" ON)
if(BUILD_CL_TOOL)
    add_executable(mscomplex-tri-cl-tool  $<TARGET_OBJECTS:mscomplex-tri-core> main.cpp)
    target_link_libraries(mscomplex-tri-cl-tool  ${Boost_LIBRARIES})
    install(TARGETS mscomplex-tri-cl-tool DESTINATION ${MSCOMPLEX_TRI_INSTALL_DIR_BIN})
endif(BUILD_CL_TOOL)

