project(mscomplex-tri-pymstri)

find_package(Eigen3 REQUIRED)

find_package(Boost 1.48 COMPONENTS python serialization system REQUIRED)

find_package(PythonLibs REQUIRED)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../core/
  ${Boost_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
)

add_library(pymstri SHARED $<TARGET_OBJECTS:mscomplex-tri-core> pymstri.cpp)

target_link_libraries(pymstri ${Boost_LIBRARIES})

set_target_properties(pymstri PROPERTIES PREFIX "")

set(PYTHON_SITE_PACKAGES_INSTALL_DIR "" CACHE PATH "installation dir for the python interface module")

if(PYTHON_SITE_PACKAGES_INSTALL_DIR)
install(TARGETS pymstri DESTINATION ${PYTHON_SITE_PACKAGES_INSTALL_DIR})
else(PYTHON_SITE_PACKAGES_INSTALL_DIR)
install(TARGETS pymstri DESTINATION ${MSCOMPLEX_TRI_INSTALL_DIR_LIB})
endif(PYTHON_SITE_PACKAGES_INSTALL_DIR)

option(INSTALL_PYMSTRI_EXAMPLES "Install the sample scripts" OFF)

if(INSTALL_PYMSTRI_EXAMPLES)
install(DIRECTORY examples/ DESTINATION 
	${MSCOMPLEX_TRI_INSTALL_DIR_SHARE})
endif(INSTALL_PYMSTRI_EXAMPLES)