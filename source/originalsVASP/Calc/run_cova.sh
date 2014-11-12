#!/bin/bash
#PBS -N cova
#PBS -j oe
#PBS -l nodes=1:ppn=10a:E5-2
#PBS -q batch
cd $PBS_O_WORKDIR
export nCores=10
pwd
vasp_alt_prim.py
vaspy
sh Calc/volume_ion.sh
sh Calc/volume_ion.sh
sh Calc/soc.sh
