#* *************************************************************#
# *  
# *   Active Particles on Curved Spaces (SAMoS)
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


if (DOXYGEN)

# policy change: always build all documentation even on platforms that don't
# have CUDA/Python/FFT/whatever
#set(DOXYGEN_PREDEFINED ${DOXYGEN_PREDEFINED} "  \\
#   \"BOOST_VERSION=103500\" ")

# disable DOT in doxygen if dot is not found
if(DOXYGEN_DOT_PATH)
	set(DOXYGEN_HAVE_DOT "YES")
else(DOXYGEN_DOT_PATH)
	set(DOXYGEN_HAVE_DOT "NO")
endif(DOXYGEN_DOT_PATH)


# configure the doxygen file
configure_file (${SAMoS_SOURCE_DIR}/doc/Doxyfile.in ${SAMoS_BINARY_DIR}/doc/Doxyfile)

add_custom_command (OUTPUT ${SAMoS_BINARY_DIR}/doc/html/index.html
					COMMAND ${DOXYGEN} ${SAMoS_BINARY_DIR}/doc/Doxyfile
					DEPENDS ${SAMoS_SRCS}
					${SAMoS_HDRS}
					${SAMoS_SOURCE_DIR}/doc/Doxyfile.in)
					

add_custom_target (doc ALL DEPENDS ${SAMoS_BINARY_DIR}/doc/html/index.html)

endif (DOXYGEN)
