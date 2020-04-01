# * *************************************************************
# *  
# *   Soft Active Mater on Surfaces (SAMoS)
# *   
# *   Author: Rastko Sknepnek
# *  
# *   Division of Physics
# *   School of Engineering, Physics and Mathematics
# *   University of Dundee
# *   
# *   (c) 2013, 2014
# * 
# *   School of Science and Engineering
# *   School of Life Sciences 
# *   University of Dundee
# * 
# *   (c) 2015
# * 
# *   Author: Silke Henkes
# * 
# *   Department of Physics 
# *   Institute for Complex Systems and Mathematical Biology
# *   University of Aberdeen  
# * 
# *   (c) 2014, 2015
# *  
# *   This program cannot be used, copied, or modified without
# *   explicit written permission of the authors.
# * 
# * ************************************************************** 


CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
if(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

project (SAMoS C CXX)


if(ENABLE_STATIC)
  SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
  SET (THREAD_LIB -lpthread)
endif(ENABLE_STATIC)

# Setup a number of misc options and libraries
include (CMakeMiscSetup.txt)
# Find the boost libraries and set them up
include (CMakeBoostSetup.txt)
# Find GSL libraries
include (CMakeGSLSetup.txt)
# Set default CFlags
include (CMakeCFlagsSetup.txt)
# Configure some source files, include directories, and create variables listing all source files
include (CMakeSRCSetup.txt)
# Configure VTK libraries
include (CMakeVTKSetup.txt)
# Configure CGAL libraries
include (CMakeCGALSetup.txt)

################################
## Define common libraries used by every target in MEMBRANE
set(BOOST_LIBS 	${Boost_REGEX_LIBRARY} ${Boost_IOSTREAMS_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} )
set(GSL_LIBS ${GSL_LIBRARIES})
set(MATH_LIB -lm)
set(SAMoS_LIBS ${GSL_LIBS} ${MATH_LIB} ${BOOST_LIBS} ${VTK_LIBS})



# ##############################################
# place all executables in the build directory 
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

###############################
## include documentation directories
if (ENABLE_DOXYGEN)
	add_subdirectory (doc)
endif (ENABLE_DOXYGEN)

add_subdirectory(src)


