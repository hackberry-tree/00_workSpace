#!/bin/bash
#PBS -N phonon_job
#PBS -j oe
#PBS -l nodes=1:ppn=16
#PBS -q default
cd $PBS_O_WORKDIR

/home/ohtani/local/openmpi-1.8.1-intel64/bin/mpirun -n 16 /home/ohtani/local/vasp.5.2/vasp
