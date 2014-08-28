#!/bin/bash
cp KPOINTS_relax_reduced KPOINTS
cp INCAR_cell INCAR

mpirun -n $nCores /opt/vasp5/vasp.5.2/vasp
judgeRX.py

clearFiles.sh