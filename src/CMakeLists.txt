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
target_link_libraries(samos  ${SAMoS_LIBS} ${THREAD_LIB})


set_target_properties(samos PROPERTIES PREFIX ""  OUTPUT_NAME "samos")

