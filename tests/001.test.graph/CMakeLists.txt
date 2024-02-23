cmake_minimum_required(VERSION 3.13.4)
set( CMAKE_VERBOSE_MAKEFILE ON )
#you always have to create a project. This is what your files
#are named, how you create objects, etc.
project(glygraph VERSION 0.0.1 DESCRIPTION "GLYCAM Graph Library")

################## BEGIN GETTING TOOLING TOGETHER ###################

#Exit if user doesnt have a new enough compiler to use c++17
#Need to check its at least version 7.0 cause that is supposed to
#be when c++17 was full supported by g++
if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 7.0)
    message(WARNING "Compiler used: ${CMAKE_CXX_COMPILER_ID}")
    message(WARNING "Compiler version used: ${CMAKE_CXX_COMPILER_VERSION}")
    message(WARNING "Error, your compiler version is a high enough version")
    #fatal error msg causes it to stop and exit
    message(FATAL_ERROR "Check gcc version is at least 7.0")
endif()

#Show user what compiler they using
message("################## TOOLS CMAKE IS USING ###################\n")
message(STATUS "CMake Version:\t${CMAKE_VERSION}")
message("")
#https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html
message(STATUS "C++ Compiler Used:\t${CMAKE_CXX_COMPILER_ID}")
message(STATUS "C++ Compiler Ver:\t${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "C++ Compiler Path:\t${CMAKE_CXX_COMPILER}")
message("")
message(STATUS "C Compiler Used:\t${CMAKE_C_COMPILER_ID}")
message(STATUS "C Compiler Version:\t${CMAKE_C_COMPILER_VERSION}")
message(STATUS "C Compiler Path:\t${CMAKE_C_COMPILER}")
message("")

################## END GETTING TOOLING TOGETHER ###################

################## BEGIN SETTING BASE COMPILE FLAGS ###################

#export compile commands for sourcetrail, etc.
#it is just a json of how commands interact, most tools use this
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#Now we want to override the typical build type flags. Original ones are as follows:
#1. Release: `-O3 -DNDEBUG`
#2. Debug: `-O0 -g`
#3. RelWithDebInfo: `-O2 -g -DNDEBUG`
#4. MinSizeRel: `-Os -DNDEBUG`
#NOTE: Recently when messing with some other
#cmake stuff, found out that these do not necessarily
#override the default commands. Need to fix this....
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O0 -DNDEBUG")

################## END SETTING BASE COMPILE FLAGS ###################
#add a library that will be an object libary built out of all our cfiles. Dunno if this should stay shared
add_library(${PROJECT_NAME} SHARED graph.cpp)

#help to have abs paths
get_filename_component(GMML_BASE_DIR "${CMAKE_SOURCE_DIR}../../../" ABSOLUTE)

#Set up paths to gmml.so lib stuff
add_library(gmml SHARED IMPORTED PUBLIC)
set_target_properties(gmml PROPERTIES IMPORTED_LOCATION "${GMML_BASE_DIR}/lib/libgmml.so")

#set up our header sources for the executable
set(HEADER_FILES
    ${GMML_BASE_DIR}/includes/MolecularModeling/TemplateGraph/Algorithms/include/SubgraphMatching.hpp
    ${GMML_BASE_DIR}/includes/MolecularModeling/TemplateGraph/Algorithms/include/TotalCycleDecomposition.hpp
    ${GMML_BASE_DIR}/includes/MolecularModeling/TemplateGraph/GraphStructure/include/Edge.hpp
    ${GMML_BASE_DIR}/includes/MolecularModeling/TemplateGraph/GraphStructure/include/Graph.hpp
    ${GMML_BASE_DIR}/includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp
    ${GMML_BASE_DIR}/includes/MolecularModeling/TemplateGraph/LazyPrints/LazyPrinters.hpp
    ${GMML_BASE_DIR}/includes/MolecularModeling/TemplateGraph/Algorithms/include/ConnectivityIdentifier.hpp
    )


add_executable(test graph.cpp ${HEADER_FILES})
target_link_libraries(test gmml)

include_directories("${GMML_BASE_DIR}")

target_compile_options(test PUBLIC -Wall -W)

#multithread supp
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

#Actually have our testing
#enable_testing()
#add_test(NAME test_test
#            COMMAND cmp test "${CMAKE_SOURCE_DIR}/graphTest.txt")
