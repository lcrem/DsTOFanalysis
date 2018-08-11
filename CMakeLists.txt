# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(dstof)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Net)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
ROOT_GENERATE_DICTIONARY(G__RawDsTofHeader RawDsTofHeader.h LINKDEF LinkDef.h)

#---Create a shared library with geneated dictionary
add_library(RawDsTofHeader SHARED RawDsTofHeader.cxx G__RawDsTofHeader.cxx)
target_link_libraries(RawDsTofHeader ${ROOT_LIBRARIES})

#---Create  a main program using the library
add_executable(makeNiceTree makeNiceTree.cxx)
target_link_libraries(makeNiceTree RawDsTofHeader)
add_executable(makeHitMap makeHitMap.cxx)
target_link_libraries(makeHitMap RawDsTofHeader)
