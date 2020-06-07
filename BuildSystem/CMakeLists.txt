cmake_minimum_required(VERSION 3.12)
project(IBM)

# Variables
SET(EXE_NAME "ibm")

# Change the output folders for our build to something cleaner and platform dependent
SET(CMAKE_ROOT   bin/cmake)

STRING(REPLACE "BuildSystem" "Source" PROJECT_HOME_DIRECTORY ${CMAKE_HOME_DIRECTORY})

message("Trace for OpenMP_CXX_FLAGS variable: ${OpenMP_CXX_FLAGS}")

# GLOBAL BUILD VARIABLES AND LIBRARIES
IF (DEBUG)
 SET(COMPILE_OPTIONS " -O0 -g2 -Wall -Werror -fmessage-length=0 -Wno-unused-local-typedefs")
ELSEIF (RELEASE)
 SET(COMPILE_OPTIONS " -O2 -g0 -Wall -Werror -fmessage-length=0 -Wno-unused-local-typedefs")
ELSE()
 SET(COMPILE_OPTIONS " -O2 -g0 -Wall -Werror -fmessage-length=0 -Wno-unused-local-typedefs")
ENDIF ()

OPTION (USE_OpenMP "Use OpenMP" ON)
IF(USE_OpenMP)
  FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
    set (COMPILE_OPTIONS "${COMPILE_OPTIONS} ${OpenMP_CXX_FLAGS}")
  ELSE()
    message(FATAL_ERROR "Could not find OpenMP on your compiler, please check that you have a compiler that supports OpenMP")
  ENDIF()
ENDIF()

# Set variables based on OS
IF (WIN32)
 ## WINDOWS BUILD VARIABLES AND LIBRARIES
 SET(buildOS windows) 
 SET(COMPILE_OPTIONS "${COMPILE_OPTIONS} -D__GNUDOS__ -std=c++14")
 SET(thirdPartyLibraries ${thirdPartyLibraries} pthread)
 
ELSE (WIN32)
 ## LINUX BUILD VARIABLES AND LIBRARIES
 SET(buildOS linux)
 SET(COMPILE_OPTIONS "${COMPILE_OPTIONS} -fPIC -std=c++14")
 SET(thirdPartyLibraries ${thirdPartyLibraries} pthread rt)
ENDIF ()

  
# BUILD PARAMETERS

MESSAGE("COMPILE_OPTIONS: ${COMPILE_OPTIONS}")
MESSAGE("SPECIAL_LIBRARIES: ${specialLibraries}")

INCLUDE_DIRECTORIES(SYSTEM "${CMAKE_HOME_DIRECTORY}/bin/${buildOS}/thirdparty/include")
INCLUDE_DIRECTORIES("${PROJECT_HOME_DIRECTORY}")

# This snippet of code will find all of our source and test files
# and auto-populate variables with them for the build
FILE(GLOB_RECURSE sourceFiles ${PROJECT_HOME_DIRECTORY}/*.cpp)


# Get a list of our third party libraries
IF (RELEASE)
 SET(thirdPartyFolder "bin/${buildOS}/thirdparty/lib/release/")
ELSE (RELEASE)
 SET(thirdPartyFolder "bin/${buildOS}/thirdparty/lib/debug/")
ENDIF ()
MESSAGE("Third Party Folder: ${thirdPartyFolder}")
LINK_DIRECTORIES("${CMAKE_HOME_DIRECTORY}/${thirdPartyFolder}")
LINK_DIRECTORIES("${CMAKE_HOME_DIRECTORY}/bin/${buildOS}/thirdparty/lib/special/")

FILE(GLOB thirdPartyLibrariesRaw RELATIVE ${CMAKE_HOME_DIRECTORY} ${thirdPartyFolder}*.a)
FOREACH(library ${thirdPartyLibrariesRaw})
  STRING(REPLACE ${thirdPartyFolder} "" shortLibrary ${library})
  SET(thirdPartyLibraries ${thirdPartyLibraries} -l:${shortLibrary})
ENDFOREACH()

FOREACH(library ${specialLibraries})
SET(thirdPartyLibraries ${thirdPartyLibraries} ${library})
ENDFOREACH()

SET(LINK_OPTIONS "${OpenMP_CXX_FLAGS}")

IF (LIBRARY)
  ADD_LIBRARY(${EXE_NAME} SHARED ${sourceFiles} )
  IF(WIN32)
    SET(LINK_OPTIONS "-static")
  ENDIF(WIN32)
ELSE(LIBRARY)
  MESSAGE("Link Options: ${LINK_OPTIONS}")
  MESSAGE("Third party libraries: ${thirdPartyLibraries}")

  ADD_EXECUTABLE(${EXE_NAME} ${sourceFiles} )
ENDIF(LIBRARY)

SET_TARGET_PROPERTIES(${EXE_NAME} PROPERTIES COMPILE_FLAGS ${COMPILE_OPTIONS} LINK_FLAGS ${LINK_OPTIONS})
TARGET_LINK_LIBRARIES(${EXE_NAME} ${thirdPartyLibraries})
