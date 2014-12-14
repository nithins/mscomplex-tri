project(mscomplex-tri-cli)

set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(Boost COMPONENTS program_options serialization system REQUIRED)

find_package(Eigen3 REQUIRED)

set(SRCS
  trimesh.h

  trimesh_dataset.h
  trimesh_dataset_ensure.h
  trimesh_dataset.cpp

  trimesh_mscomplex.h
  trimesh_mscomplex_ensure.h
  trimesh_mscomplex_simp.h
  trimesh_mscomplex.cpp

  tri_edge.h
  tri_edge.cpp

  utl.h
  utl.cpp

  main.cpp
  )


include_directories(
  ${EIGEN3_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${Boost_INCLUDE_DIRS})

add_executable(${PROJECT_NAME}  ${SRCS})

target_link_libraries(${PROJECT_NAME}  ${Boost_LIBRARIES})


