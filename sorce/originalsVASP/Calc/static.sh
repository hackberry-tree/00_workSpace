#!/bin/bash
cp KPOINTS_relax KPOINTS
cp INCAR_static INCAR

mpirun -n $nCores /opt/vasp5/vasp.5.2/vasp
cp OSZICAR osz.static
cp OUTCAR out.static
cp vasprun.xml xml.static

clearFiles.sh
