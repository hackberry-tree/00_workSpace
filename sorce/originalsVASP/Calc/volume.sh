#!/bin/bash
cp KPOINTS_relax KPOINTS
cp INCAR_volume INCAR

mpirun -n $nCores /opt/vasp5/vasp.5.2/vasp
judgeRX.py

clearFiles.sh