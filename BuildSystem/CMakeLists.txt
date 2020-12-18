cmake_minimum_required(VERSION 3.12)
project(IBM)

# Variables
SET(EXE_NAME "ibm")

# Change the output folders for our build to something cleaner and platform dependent
SET(CMAKE_ROOT   bin/cmake)

STRING(REPLACE "BuildSystem" "Source" PROJECT_HOME_DIRECTORY ${CMAKE_HOME_DIRECTORY})

# GLOBAL BUILD VARIABLES AND LIBRARIES
IF (DEBUG)
 SET(COMPILE_OPTIONS " -O0 -g2 -Wall -Werror -fmessage-length=0 -Wno-unused-local-typedefs -fopenmp")
ELSEIF (RELEASE)
 SET(COMPILE_OPTIONS " -O2 -g0 -Wall -Werror -fmessage-length=0 -Wno-unused-local-typedefs -fopenmp")
ELSE()
 SET(COMPILE_OPTIONS " -O2 -g0 -Wall -Werror -fmessage-length=0 -Wno-unused-local-typedefs -fopenmp")
ENDIF ()
########################################################################################
# Rinside calls, now requires R.exe to be in path
SET(NUM_TRUNC_CHARS 2)
## Set R cppflags variables
message( "Starting RInside cmake component, if this fails install the R-library install.packages(c('RInside','Rcpp'))")
EXECUTE_PROCESS(COMMAND R CMD config --cppflags OUTPUT_VARIABLE RCPPFLAGS)

EXECUTE_PROCESS(COMMAND R CMD config --cppflags OUTPUT_VARIABLE RCPPFLAGS)
message( "RCPPFLAGS ${RCPPFLAGS}" )
## Remove -I 
STRING(REGEX REPLACE "-I" "" RCPPFLAGS ${RCPPFLAGS})
#message( "RCPPFLAGS ${RCPPFLAGS}" )
string(REPLACE " " ";" RCPPFLAGS_LIST ${RCPPFLAGS})
#message( "RCPPFLAGS_LIST ${RCPPFLAGS_LIST}" )
include_directories(${RCPPFLAGS_LIST})
#include_directories(${RCPPFLAGS})

## ldflags variables
EXECUTE_PROCESS(COMMAND R CMD config --ldflags OUTPUT_VARIABLE RLDFLAGS)
#message( "RLDFLAGS ${RLDFLAGS}" )
## get library path
if (${RLDFLAGS} MATCHES "[-][L]([^ ;])+")
    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RLDFLAGS_L)
    string(STRIP ${RLDFLAGS_L} RLDFLAGS_L )
    link_directories(${RLDFLAGS_L})
endif()
#message( "RLDFLAGS_L ${RLDFLAGS_L}" )
## get actual libraries
if (${RLDFLAGS} MATCHES "[-][l]([^;])+")
    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RLDFLAGS_l)
    STRING(REGEX REPLACE "-l" "" RLDFLAGS_l ${RLDFLAGS_l})
    string(STRIP ${RLDFLAGS_l} RLDFLAGS_l )
endif()
string(REPLACE " " ";" RLDFLAGS_l_LIST ${RLDFLAGS_l})
#message( "RCPPFLAGS_LIST ${RLDFLAGS_l_LIST}" )
#message( "RLDFLAGS_l ${RLDFLAGS_l}" )

##Rcpp
execute_process(COMMAND Rscript -e "Rcpp:::CxxFlags()" OUTPUT_VARIABLE RCPPINCL)
string(SUBSTRING ${RCPPINCL} ${NUM_TRUNC_CHARS} -1 RCPPINCL)
include_directories(${RCPPINCL})
#message( "RCPPINCL ${RCPPINCL}" )
execute_process(COMMAND Rscript -e "Rcpp:::LdFlags()" OUTPUT_VARIABLE RCPPLIBS)

if (${RCPPLIBS} MATCHES "[-][L]([^ ;])+")
    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RCPPLIBS_L)
    link_directories(${RCPPLIBS_L} )
endif()

if (${RCPPLIBS} MATCHES "[-][l][R]([^;])+")
    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RCPPLIBS_l)
endif()
#message( "RCPPLIBS ${RCPPLIBS}" )

## RInside 
execute_process(COMMAND Rscript -e "RInside:::CxxFlags()"
                OUTPUT_VARIABLE RINSIDEINCL)
string(SUBSTRING ${RINSIDEINCL} ${NUM_TRUNC_CHARS} -1 RINSIDEINCL)
include_directories(${RINSIDEINCL})
#message( "RINSIDEINCL ${RINSIDEINCL}" )

STRING(REGEX REPLACE "include" "lib/x64" RINSIDE_LIB ${RINSIDEINCL})
#message( "RINSIDE_LIB ${RINSIDE_LIB}" )
link_directories(${RINSIDE_LIB})


execute_process(COMMAND Rscript -e "RInside:::LdFlags()" OUTPUT_VARIABLE RINSIDELIBS)
if (${RINSIDELIBS} MATCHES "[-][L]([^ ;])+")
    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RINSIDELIBS_L)
    link_directories(${RINSIDELIBS_L})
endif()
#message( "RINSIDELIBS_L ${RINSIDELIBS_L}" )

if (${RINSIDELIBS} MATCHES "[-][l]([^;])+")
    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RINSIDELIBS_l)
    STRING(REGEX REPLACE "-l" "" RINSIDELIBS_l ${RINSIDELIBS_l})
endif()
#message( "RINSIDELIBS_l ${RINSIDELIBS_l}" )
#message( "RINSIDELIBS ${RINSIDELIBS}" )

## BLAS
execute_process(COMMAND R CMD config BLAS_LIBS OUTPUT_VARIABLE RBLAS)
## get library path
if (${RBLAS} MATCHES "[-][L]([^ ;])+")
    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RBLAS_L)
    string(STRIP ${RBLAS_L} RBLAS_L )
    ##link_directories(${RLDFLAGS_L})## same link directory as RLDFLAGS
endif()
#message( "RBLAS_L ${RBLAS_L}" )
## get actual libraries
if (${RBLAS} MATCHES "[-][l]([^;])+")
    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RBLAS_l)
    STRING(REGEX REPLACE "-l" "" RBLAS_l ${RBLAS_l})
    string(STRIP ${RBLAS_l} RBLAS_l )
endif()
#message( "RBLAS_l ${RBLAS_l}" )

##LAPACK
execute_process(COMMAND R CMD config LAPACK_LIBS OUTPUT_VARIABLE RLAPACK)
## get library path
if (${RLAPACK} MATCHES "[-][L]([^ ;])+")
    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RLAPACK_L)
    string(STRIP ${RLAPACK_L} RLAPACK_L )
    ##link_directories(${RLAPACK_L})## same link directory as RLDFLAGS
endif()
#message( "RLAPACK_L ${RLAPACK_L}" )
## get actual libraries
if (${RLAPACK} MATCHES "[-][l]([^;])+")
    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RLAPACK_l)
    STRING(REGEX REPLACE "-l" "" RLAPACK_l ${RLAPACK_l})
    string(STRIP ${RLAPACK_l} RLAPACK_l )
endif()
#message( "RLAPACK_l ${RLAPACK_l}" )

## set cmake cxx flags 
execute_process(COMMAND R CMD config CXXFLAGS OUTPUT_VARIABLE RCXXFLAGS)
#message( "RCXXFLAGS ${RCXXFLAGS}" )

#message( "CMAKE_CXX_FLAGS before changing ${CMAKE_CXX_FLAGS}" )
#set(CMAKE_CXX_FLAGS "-s -Wall ${CMAKE_CXX_FLAGS} ${RCXXFLAGS}")
#add_definitions(${CMAKE_CXX_FLAGS})

##SET(COMPILE_OPTIONS "${COMPILE_OPTIONS} ")




########################################################################################
# Set variables based on OS
IF (WIN32)
 ## WINDOWS BUILD VARIABLES AND LIBRARIES
 SET(buildOS windows) 
 SET(COMPILE_OPTIONS "${COMPILE_OPTIONS} -D__GNUDOS__ -std=c++14 -s ${RCXXFLAGS}")
 SET(thirdPartyLibraries ${thirdPartyLibraries} pthread)
 
ELSE (WIN32)
 ## LINUX BUILD VARIABLES AND LIBRARIES
 SET(buildOS linux)
 SET(COMPILE_OPTIONS "${COMPILE_OPTIONS} -std=c++14")
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

SET(LINK_OPTIONS "-fopenmp")

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
TARGET_LINK_LIBRARIES(${EXE_NAME} ${RLAPACK_l})
TARGET_LINK_LIBRARIES(${EXE_NAME} ${RBLAS_l})
TARGET_LINK_LIBRARIES(${EXE_NAME} ${RCPPLIBS_l})
TARGET_LINK_LIBRARIES(${EXE_NAME} "RInside")
TARGET_LINK_LIBRARIES(${EXE_NAME} ${RLDFLAGS_l_LIST})

## Print all the libraries
message( "Link libraries -l" )
get_target_property(OUT ${EXE_NAME} LINK_LIBRARIES)
foreach (_variableName ${OUT})
    message(STATUS "lib='${_variableName}'")
endforeach()
## Print all the library directories
message( "Link directories -L" )
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY LINK_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "dir='${dir}'")
endforeach()

