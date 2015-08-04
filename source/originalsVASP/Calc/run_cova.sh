#!/bin/bash
#PBS -N cova
#PBS -j oe
#PBS -l nodes=1:ppn=10a:E5-2
#PBS -q batch
cd $PBS_O_WORKDIR
export nCores=10
pwd
vasp_poscar.py --prim
cp POSCAR POSCAR_orig
cp POSCAR_prim POSCAR
vaspy
sh Calc/volume_ion.sh
sh Calc/volume_ion.sh
sh Calc/soc.sh
