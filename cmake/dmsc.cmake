find_package(Boost COMPONENTS program_options thread REQUIRED)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(Eigen3 REQUIRED)

set(SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh.h

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh_dataset.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh_dataset_ensure.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh_dataset.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh_mscomplex.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh_mscomplex_ensure.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh_mscomplex.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/utls/include/cpputils.h
  ${CMAKE_CURRENT_SOURCE_DIR}/utls/src/cpputils.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/utls/include/tri_edge.h
  ${CMAKE_CURRENT_SOURCE_DIR}/utls/src/tri_edge.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/utls/include/timer.h
  ${CMAKE_CURRENT_SOURCE_DIR}/utls/src/timer.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/main.cpp
  )


include_directories(
  ${EIGEN3_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/
  ${CMAKE_CURRENT_SOURCE_DIR}/utls/include/
  ${Boost_INCLUDE_DIRS})

add_executable(${PROJECT_NAME}  ${SRCS})

target_link_libraries(${PROJECT_NAME}  ${Boost_LIBRARIES})


