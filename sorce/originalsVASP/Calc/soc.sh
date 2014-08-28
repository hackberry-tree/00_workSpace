#!/bin/bash

#non-collinear
cp INCAR_presoc_nc INCAR
cp KPOINTS_soc KPOINTS
mpirun -n $nCores /opt/vasp5n/vasp.5.2/vasp
cp OSZICAR OSZICAR_presoc_nc
cp OUTCAR OUTCAR_presoc_nc

#soc
cp INCAR_soc001 INCAR
mpirun -n $nCores /opt/vasp5n/vasp.5.2/vasp
cp OUTCAR OUTCAR_soc001
cp OSZICAR OSZICAR_soc001

cp INCAR_soc100 INCAR
mpirun -n $nCores /opt/vasp5n/vasp.5.2/vasp
cp OUTCAR OUTCAR_soc100
cp OSZICAR OSZICAR_soc100

cp INCAR_soc110 INCAR
mpirun -n $nCores /opt/vasp5n/vasp.5.2/vasp
cp OUTCAR OUTCAR_soc110
cp OSZICAR OSZICAR_soc110

cp INCAR_soc111 INCAR
mpirun -n $nCores /opt/vasp5n/vasp.5.2/vasp
cp OUTCAR OUTCAR_soc111
cp OSZICAR OSZICAR_soc111

clearFiles.sh
