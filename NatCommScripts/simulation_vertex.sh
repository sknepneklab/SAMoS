#!/bin/bash
# *****************************************************************************
# *
# *  This BASH script is a part of tha analysis of the data published in 
# *  the paper: "Universal motion patterns in confluent cell monolayers"
# *  by Silke Henkes, Kaja Kostanjevec, J. Martin Collinson, Rastko Sknepnek, 
# *  and Eric Bertin, Jounral name, vol, page (2019).
# *
# *  Please refer to the document Computational_summary.pdf for a detailed
# *  description of the tasks performed by this script.
# * 
# *****************************************************************************


if ! [ -x "$(command -v samos)" ]; then
  echo 'Error: samos exacutable cannot be found.' >&2
  echo '       Please make sure to add the path to the samos executable' >&2
  echo '       to the PATH shell variable.' >&2
  exit 1
fi

topdir_code=`command -v samos`

if [ -z "$OUTDIR" ]; then
  echo "Warning!" >&2
  echo "Output directory has not been set. Using current directory." >&2
  echo "You can set output directory by setting the shell environment variable OUTDIR." >&2
  echo "In bash, e.g., you can type export OUTDIR=/path/to/output" s>&2
  topdir_out="./"
else
  topdir_out=$OUTDIR
fi

topdir=$topdir_out
parfile='makelarge_vertex.conf'
equilfile='equilibration_vertex.conf'

# Here, vary the active temperature v0^2/(nu_r), at fixed temperature and then vary tau
tempval='0.1 0.05 0.02 0.01 0.005 0.002'
# active noise level, \tau = 2/nu
nuval='10.0 1.0 0.1 0.01 0.001'
# demonstration values
tempval='0.01'
nuval='0.01'

# runtime in units of time steps
trun=5000000
# short demonstration run
trun=200000

# for example
seed=367496
for temp in ${tempval}
do
	echo ${temp}
	for nu in ${nuval}
	do
    echo "nu: " ${nu}
		# convert the active temperature value into a driving velocity, T_eff = v_0^2 /nu
		v0=`echo "sqrt(${temp}*${nu})" | bc -l`
		echo "v0: " ${v0}
		cd ${topdir}
		confdir=${topdir_out}/vertex/data_temp_${temp}/data_nu_${nu}_plane/
		mkdir -p ${confdir}
		# create input file for the soft particle equilibration which serves as samos vertex input file
		# exact same number as soft particle system, 30% polydispersity at packing fraction 1
		# target area per cell will be pi
		python plane_circle_epithelial.py -n 3183 -p 0.3 -f 1.0 -a 3.141592 -o  ${confdir}epithelial_randomini.dat -i ${confdir}vertex_bound.dat -v 0.0 -b 0 -d 1.0
		
		# equilibrate inside the configuratin directory
		equilparfile=${confdir}${equilfile}
		cp ${equilfile} ${equilparfile}
		cd ${confdir}
		${topdir_code}samos ${equilfile}
                
		# use last equilibrated values as input for the vertex simulation
		mv equil_0000100000.dat vertex_input.dat
 		
 		# copy configuration file into run directory and edit parameters using sed
		newparfile=vertex_temp_${temp}_nu_${nu}.conf
		cp ${topdir}${parfile} ${newparfile}
		#pair_potential vp { K = 1.0; gamma = 1.0; lambda = -7.40 }
		sed -i "s|@KVAL|1.0|" $newparfile
		sed -i "s|@GAMMA|1.0|" $newparfile
		# p0 = 3.6, essentially. Note that this is lambda = - 2 P0 gamma, except there is no 2 because of how it's coded
		# then p0 = P0 / sqrt{A0) = P0 / sqrt(pi), so finally, lambda = - gamma * p0 * sqrt(pi)
		sed -i "s|@LAMBDA|-6.2|" $newparfile 
		# line tension to avoid boundary instability
		sed -i "s|@LINE|0.3|" $newparfile 
		sed -i "s|@V0|$v0|" $newparfile
		sed -i "s|@NU|$nu|" $newparfile
		sed -i "s|@TRUN|$trun|" $newparfile
		seed=$((seed+1))
		sed -i "s|@SEED|$seed|" $newparfile
		echo ${topdir_in}samos ${newparfile}
		${topdir_code}samos ${newparfile}
		
	done
done
