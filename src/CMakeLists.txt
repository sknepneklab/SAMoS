#* *************************************************************#
# *  
# *   Active Particles on Curved Spaces (APCS)
# *   
# *   Author: Rastko Sknepnek
# *  
# *   Division of Physics
# *   School of Engineering, Physics and Mathematics
# *   University of Dundee
# *   
# *   (c) 2013
# *   
# *   This program cannot be used, copied, or modified without
# *   explicit permission of the author.
# * 
# * ************************************************************* */

macro(add_source_dir dir)
add_subdirectory(${dir})

# find all source files in this directory
file(GLOB APCS_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.cpp)
file(GLOB APCS_HDRS ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.hpp)

# add them to the source group
source_group(${dir} FILES ${APCS_SRCS} ${APCS_HDRS})

# add them to the main apcs list
list(APPEND _apcs_sources ${APCS_SRCS} ${APCS_HDRS})
endmacro(add_source_dir)
#################

set(_apcs_sources "")

# process all subdirectories, which add source files to _apcs_sources 
set(_dir_list constraints  dump log  integrators  messenger  parser  potentials potentials/pair potentials/external aligner aligner/pair aligner/external system  utils)

foreach (dir ${_dir_list})
add_source_dir(${dir})
endforeach (dir)


add_executable(apcs apcs.cpp ${_apcs_sources})
target_link_libraries(apcs  ${APCS_LIBS} ${THREAD_LIB})


set_target_properties(apcs PROPERTIES PREFIX ""  OUTPUT_NAME "apcs")

