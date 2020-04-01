# * ***************************************************************************
# *
# *  Copyright (C) 2013-2016 University of Dundee
# *  All rights reserved. 
# *
# *  This file is part of SAMoS (Soft Active Matter on Surfaces) program.
# *
# *  SAMoS is free software; you can redistribute it and/or modify
# *  it under the terms of the GNU General Public License as published by
# *  the Free Software Foundation; either version 2 of the License, or
# *  (at your option) any later version.
# *
# *  SAMoS is distributed in the hope that it will be useful,
# *  but WITHOUT ANY WARRANTY; without even the implied warranty of
# *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# *  GNU General Public License for more details.
# *
# *  You should have received a copy of the GNU General Public License
# *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
# *
# * ****************************************************************************

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

target_compile_options(samos PUBLIC $<$<CONFIG:RELEASE>:-O3 -funroll-loops -ffast-math -DNDEBUG>)
target_compile_options(samos PUBLIC $<$<CONFIG:DEBUG>:-O0 -g3 -Wall>)

target_link_libraries(samos  ${SAMoS_LIBS} ${THREAD_LIB})


set_target_properties(samos PROPERTIES 
    PREFIX ""  
    OUTPUT_NAME "samos" 
    CXX_STANDARD ${CXX_STANDARD}
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    )

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX $ENV{HOME}/samos CACHE PATH "Setting default install path" FORCE)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

file(GLOB ANALYSIS_FILES "${PROJECT_SOURCE_DIR}/utils/*.py")
file(GLOB EXAMPLE1_FILES "${PROJECT_SOURCE_DIR}/configurations/cells_fixed/*")
file(GLOB EXAMPLE2_FILES "${PROJECT_SOURCE_DIR}/configurations/circular_patch/*")
file(GLOB EXAMPLE3_FILES "${PROJECT_SOURCE_DIR}/configurations/growth/*")

install(TARGETS samos DESTINATION bin)
install(FILES ${ANALYSIS_FILES} DESTINATION "analysis")
install(FILES ${EXAMPLE1_FILES} DESTINATION "examples/cells_fixed")
install(FILES ${EXAMPLE2_FILES} DESTINATION "examples/circular_patch")
install(FILES ${EXAMPLE3_FILES} DESTINATION "examples/growth")
install(DIRECTORY ${PROJECT_SOURCE_DIR}/doc/tutorial/ DESTINATION "doc/tutorial")
install(CODE "MESSAGE(\"
-------------------------------------------------------------------------------------------- 

 Please add ${CMAKE_INSTALL_PREFIX}/bin to the PATH variable 
        and ${CMAKE_INSTALL_PREFIX}/analysis to the PYTHONPATH variable 
      
 In BASH you can do this by typing: 

        export PATH=$PATH:${CMAKE_INSTALL_PREFIX}/bin \n
        export PYTHONPATH=$PYTHONPATH:${CMAKE_INSTALL_PREFIX}/analysis 
        
 In order to make the change permanent, please add above two lines to you .bashrc script. 

-------------------------------------------------------------------------------------------- \")" )