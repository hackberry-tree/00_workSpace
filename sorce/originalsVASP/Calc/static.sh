#!/bin/bash
cp KPOINTS_relax KPOINTS
cp INCAR_static INCAR

mpirun -n $nCores /opt/vasp5/vasp.5.2/vasp

clearFiles.sh
