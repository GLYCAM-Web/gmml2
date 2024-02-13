cmake_minimum_required(VERSION 3.13.4)

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

add_library(gmml SHARED IMPORTED PUBLIC)
set_target_properties( gmml PROPERTIES IMPORTED_LOCATION "${CMAKE_SOURCE_DIR}/../../lib/libgmml.so")

#add a library that will be an object libary built out of all our cfiles. Dunno if this should stay shared
add_library(${PROJECT_NAME} SHARED graph.cpp)

target_compile_options(${PROJECT_NAME} PUBLIC -Wall -W)

target_link_libraries(${PROJECT_NAME} gmml)
target_link_libraries(${PROJECT_NAME} stdc++fs)

#tryna figure out wtf is happening.....
get_property(HALP TARGET gmml PROPERTY LOCATION)
message(STATUS ${HALP})

add_executable(testDriver graph.cpp)
#lole need to do this proper.
#Currently thinking of having files read in as one arg then
#other arg(s) are the algo/tests to run
#enable_testing()
#add_test(NAME testDriver_test
#            COMMAND cmp testDriver "${CMAKE_SOURCE_DIR}/graphTest.txt")
