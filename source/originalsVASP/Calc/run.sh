#!/bin/bash
#PBS -N vasp
#PBS -j oe
#PBS -l nodes=1:ppn=4:W
#PBS -q default
cd $PBS_O_WORKDIR
export nCores=4

pwd

#sh Calc/band.sh

sh Calc/soc.sh

#sh Calc/polarized.sh

#sh Calc/volume.sh

#sh Calc/cell.sh

#sh Calc/static.sh
