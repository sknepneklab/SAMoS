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
  echo "Warning!"
  echo "Output directory has not been set. Using current directory." >&2
  echo "You can set output directory by setting the shell environment variable OUTDIR." >&2
  echo "In bash, e.g., you can type export OUTDIR=/path/to/output" >&2
  topdir_out="./"
else
  topdir_out=$OUTDIR
fi

topdir=$topdir_out
parfile='makelarge.conf'

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

# running script, including creation of directory, initial configuration and adapted paramter file
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
		confdir=${topdir_out}/soft/data_temp_${temp}/data_nu_${nu}_plane/
		mkdir -p ${confdir}
		# create input file of random disks with packing fractin 1 and 30% polydispersity
		python random_plane_poly.py -x 100 -y 100 -f 1.0 -v 0.0 -p 0.3 -o ${confdir}plane_temp_${temp}_nu_${nu}.txt
		# copy configuration file into run directory and edit parameters using sed
		newparfile=${confdir}plane_temp_${temp}_nu_${nu}.conf
		cp ${parfile} ${newparfile}
		sed -i "s|@INPUT|plane_temp_${temp}_nu_${nu}.txt|" $newparfile
		sed -i "s|@PREFIX|plane_temp_${temp}_nu_${nu}|" ${newparfile}
		sed -i "s|@V0|$v0|" $newparfile
		sed -i "s|@NU|$nu|" $newparfile
		sed -i "s|@RUN|$trun|" $newparfile
		cd ${confdir}
		# run 
		echo ${topdir_code}samos ${newparfile}
		${topdir_code}samos ${newparfile}
	done
done
