# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 2.8.10 FATAL_ERROR)
project(dstof)
set(libname "dstof")

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Net)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})


file(GLOB HEADER_FILES "Raw*.h" )

file(GLOB_RECURSE SOURCE_FILES Raw*.cxx)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
set(DICTNAME G__${libname})

ROOT_GENERATE_DICTIONARY(${DICTNAME} ${HEADER_FILES} LINKDEF LinkDef.h)

#---Create a shared library with geneated dictionary
add_library(${libname} SHARED ${SOURCE_FILES} ${DICTNAME}.cxx)
target_link_libraries(${libname} ${ROOT_LIBRARIES})


if( ${ROOT_VERSION} VERSION_GREATER "5.99.99")
    add_custom_target(${DICTNAME}.pcm DEPENDS ${DICTNAME})
endif()

message("DICTNAME is set to ${DICTNAME}" )
message("HEADER_FILES is set to ${HEADER_FILES}" )
message("SOURCE_FILES is set to ${SOURCE_FILES}" )

#---Create  a main program using the library
add_executable(makeNiceTree makeNiceTree.cxx) 
target_link_libraries(makeNiceTree ${ROOT_LIBRARIES} dstof)
add_executable(makeHitMap makeHitMap.cxx) 
target_link_libraries(makeHitMap ${ROOT_LIBRARIES} dstof)
add_executable(makeHVplots makeHVplots.cxx) 
target_link_libraries(makeHVplots ${ROOT_LIBRARIES} dstof)
add_executable(makeDebugTree makeDebugTree.cxx) 
target_link_libraries(makeDebugTree ${ROOT_LIBRARIES} dstof)
add_executable(makeQuickTofPlot makeQuickTofPlot.cxx) 
target_link_libraries(makeQuickTofPlot ${ROOT_LIBRARIES} dstof)
add_executable(makeSpillDB makeSpillDB.cxx) 
target_link_libraries(makeSpillDB ${ROOT_LIBRARIES} dstof)
