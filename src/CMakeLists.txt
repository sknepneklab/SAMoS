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

macro(add_source_dir dir)
add_subdirectory(${dir})

# find all source files in this directory
file(GLOB SAMoS_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.cpp)
file(GLOB SAMoS_HDRS ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.hpp)

# add them to the source group
source_group(${dir} FILES ${SAMoS_SRCS} ${SAMoS_HDRS})

# add them to the main samos list
list(APPEND _samos_sources ${SAMoS_SRCS} ${SAMoS_HDRS})
endmacro(add_source_dir)
#################

set(_samos_sources "")

# process all subdirectories, which add source files to _samos_sources 
set(_dir_list constraints  dump log  integrators  messenger  parser  potentials potentials/pair potentials/external potentials/bond potentials/angle aligner aligner/pair aligner/external system  utils population)

foreach (dir ${_dir_list})
add_source_dir(${dir})
endforeach (dir)


add_executable(samos samos.cpp ${_samos_sources})
target_link_libraries(samos  ${SAMoS_LIBS} ${THREAD_LIB})


set_target_properties(samos PROPERTIES PREFIX ""  OUTPUT_NAME "samos")

