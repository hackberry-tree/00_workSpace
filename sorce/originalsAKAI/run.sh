#!/bin/bash
#PBS -N job
#PBS -j oe
#PBS -l nodes=1:ppn=1:curie.tagen.tohoku.ac.jp
#PBS -q default
cd $PBS_O_WORKDIR
#
/home/enoki/cpa2002v009c/specx < input > output
