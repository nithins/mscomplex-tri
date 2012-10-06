find_path(PYTHONQT_INCLUDE_DIR
  NAME PythonQt.h
  PATHS /usr/include/PythonQt /usr/local/include/PythonQt
  )

find_library(PYTHONQT_LIBRARY
  NAMES PythonQt
  PATH /usr/lib /usr/local/lib
  )

if(PYTHONQT_INCLUDE_DIR AND PYTHONQT_LIBRARY)
  set(PythonQt_FOUND TRUE)
  message(STATUS "Found PythonQt Include Directory: ${PYTHONQT_INCLUDE_DIR}")
  message(STATUS "Found PythonQt Library: ${PYTHONQT_LIBRARY}")
else(PYTHONQT_INCLUDE_DIR AND PYTHONQT_LIBRARY)
  message(FATAL_ERROR "Could not find PythonQt")
  set(PythonQt_FOUND FALSE)
endif(PYTHONQT_INCLUDE_DIR AND PYTHONQT_LIBRARY)