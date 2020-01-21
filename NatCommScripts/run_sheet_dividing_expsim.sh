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

parfile='epithelial_divide_tracers_expsim.conf'

# death rate (inverse hours per particle)
rdeathval='0.01'
nuval='0.8'
kval='55'
v0val='90'
# attempted division rate (inverse hours per particle)
multval='0.05' 

lx=867.0
ly=662.0
radius=11.0
trun=200000
freq=100

seed=328574
echo $seed


for rdeath in ${rdeathval}
do
  echo "rdeath: " ${rdeath}
	for v0 in ${v0val}
	do
		echo "v0: " $v0
		for nu in ${nuval}
		do
			echo "nu: " ${nu}
			for mult in ${multval}
			do
        echo "nu: " ${nu}
        for k in ${kval}
        do
            echo "k: " ${k}
            confdir=${topdir_out}divide/rdeath_${rdeath}/nu_${nu}/v0_${v0}/ratio_${mult}/k_${k}/conf3/
            echo "confdir: " ${confdir}
            mkdir -p ${confdir}
            rd=`echo "${rdeath}/${mult}" | bc -l`
            echo "rd: " ${rd}

            python ${topdir_out}/random_plane_poly_expsim.py -x ${lx} -y ${ly} -f 1.0 -v $v0 -p 0.3 -i ${radius} -r 0.05 -o ${confdir}plane.txt
            newparfile=${confdir}${parfile}
            
            cp ${topdir_out}${parfile} ${newparfile}
            # sed patterns for the changing variables
            sed -i "s|@KVAL|$k|" $newparfile
            sed -i "s|@NUVAL|$nu|" $newparfile
            sed -i "s|@V0|$v0|" $newparfile
            sed -i "s|@DIV|$rd|" $newparfile
            sed -i "s|@DEATH|$rdeath|" $newparfile
            sed -i "s|@FREQ|$freq|" $newparfile
            sed -i "s|@TRUN|$trun|" $newparfile
            sed -i "s|@SEED|$seed|" $newparfile
            cd $confdir
            echo ${topdir_code}samos $newparfile
            ${topdir_code}samos $newparfile

            seed=$(($seed+1))
				
				done
			done
		done
	done
done
	    
