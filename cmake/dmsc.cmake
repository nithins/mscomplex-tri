find_package(Boost COMPONENTS program_options thread REQUIRED)

set(SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh.h

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh_dataset.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh_dataset_ensure.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh_dataset.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh_mscomplex.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh_mscomplex_ensure.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh_mscomplex.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh_datamanager.h
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/trimesh_datamanager.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/utls/include/cpputils.h
  ${CMAKE_CURRENT_SOURCE_DIR}/utls/src/cpputils.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/utls/include/n_vector.h
  ${CMAKE_CURRENT_SOURCE_DIR}/utls/src/n_vector.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/utls/include/tri_edge.h
  ${CMAKE_CURRENT_SOURCE_DIR}/utls/src/tri_edge.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/utls/include/timer.h
  ${CMAKE_CURRENT_SOURCE_DIR}/utls/src/timer.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/main.cpp
  )


include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/dmsc/
  ${CMAKE_CURRENT_SOURCE_DIR}/utls/include/
  ${Boost_INCLUDE_DIRS})

add_executable(${PROJECT_NAME}  ${SRCS})

target_link_libraries(${PROJECT_NAME}  ${Boost_LIBRARIES})


