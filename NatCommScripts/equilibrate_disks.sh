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


# demonstration values
temp='0.01'
nuval='0.01'

for nu in ${nuval}
do
    echo ${nu}
    confdir=${topdir_out}/soft/data_temp_${temp}/data_nu_${nu}_plane/
    newparfile=${confdir}${parfile}
    cp ${parfile} ${newparfile}
    # in actual simulations: starting configuration is the midpoint, at 2.5 million time steps
    # inputfile=${confdir}plane_temp_${temp}_nu_${nu}_0002500000.dat
    # for demonstration purposes: the last positions at step 50000
    inputfile=${confdir}plane_temp_${temp}_nu_${nu}_0000050000.dat
    echo ${confdir}
    echo ${newparfile}
    sed -i "s|@INPUT|$inputfile|" $newparfile
    # 100000 step equilibration run with purely overdamped / steepest descent dynamics (no noise, no activity, no thermal fluctuations)
    echo ${topdir_in}samos ${newparfile}
    cd ${confdir}
    ${topdir_code}samos ${newparfile}
done
