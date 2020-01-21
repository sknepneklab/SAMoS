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

equilfile='equilibration_vertex_expsim.conf'
parfile='sheet_vertex_expsim.conf'

# final values of the parameters 
nuval='0.8'
kval='55'
v0val='90'
p0val='3.6'

# spatial dimensions in microns
lx=867.0
ly=662.0
radius=11.0
# runtime: 200000 dt, at dt = 0.001 hours
trun=200000
freq=100

seed=985379
echo $seed


for p0 in ${p0val}
do
	for v0 in ${v0val}
	do
		echo "v0: " $v0
		for nu in ${nuval}
		do
			echo "nu: " ${nu}
      for k in ${kval}
      do
                        
        echo ${k}
        cd "$topdir_out"
        confdir="${topdir_out}"/nodiv_vertex/p0_${p0}/nu_${nu}/v0_${v0}/k_${k}/conf2/
        echo ${confdir}
        mkdir -p ${confdir}
        python plane_circle_epithelial_expsim.py -n 1500 -p 0.3 -f 1.0 -a 382 -o  ${confdir}epithelial_randomini.dat -i ${confdir}vertex_bound.dat -v 0.0 -b 0 -d 1.0

        equilparfile=${confdir}${equilfile}
        cp ${equilfile} ${equilparfile}
        cd ${confdir}
        "${topdir_code}"samos ${equilfile}

        mv equil_0000020000.dat vertex_input.dat

        newparfile=vertex.conf
        cp ${topdir_out}${parfile} ${newparfile}
        #pair_potential vp { K = 1.0; gamma = 1.0; lambda = -7.40 }
        # Gamma / radius^2, radius is 11 microns
        kappa=$(echo "${k}/121.0" | bc -l)
        echo ${kappa}
        sed -i "s|@KVAL|${kappa}|" $newparfile
        sed -i "s|@GAMMA|$k|" $newparfile
        # p0 = 3.6, essentially. Note that this is lambda = - 2 P0 gamma, except there is no 2 because of how it's coded
        # then p0 = P0 / sqrt{A0) = P0 / sqrt(pi), so finally, lambda = - gamma * p0 * sqrt(A0)
        # square root A0 = sqrt(382)=19.5
        lambda=$(echo "-${k}*${p0}*19.5" | bc -l)
        echo ${lambda}
        sed -i "s|@LAMBDA|${lambda}|" $newparfile 
        sed -i "s|@LINE|5|" $newparfile 
        sed -i "s|@V0|$v0|" $newparfile
        sed -i "s|@NU|$nu|" $newparfile
        sed -i "s|@TRUN|$trun|" $newparfile
        sed -i "s|@FREQ|$freq|" $newparfile
        seed=$((seed+1))
        sed -i "s|@SEED|$seed|" $newparfile
        echo "${topdir_code}"samos ${newparfile}
        "${topdir_code}"samos ${newparfile}
      done
    done
	done
done
