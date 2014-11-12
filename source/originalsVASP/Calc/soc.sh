#!/bin/bash

#vaspy --remake_from_fin

#non-collinear
cp INCAR_presoc_nc INCAR
cp KPOINTS_soc KPOINTS
custodian_static.py
#mpirun -n $nCores /opt/vasp5n/vasp.5.2/vasp
cp OSZICAR OSZICAR_presoc_nc
cp OUTCAR OUTCAR_presoc_nc

vaspy --remake_soc

#soc
cp INCAR_soc001 INCAR
custodian_static.py
#mpirun -n $nCores /opt/vasp5n/vasp.5.2/vasp
cp OUTCAR OUTCAR_soc001
cp OSZICAR OSZICAR_soc001

cp INCAR_soc100 INCAR
custodian_static.py
#mpirun -n $nCores /opt/vasp5n/vasp.5.2/vasp
cp OUTCAR OUTCAR_soc100
cp OSZICAR OSZICAR_soc100

cp INCAR_soc110 INCAR
custodian_static.py
#mpirun -n $nCores /opt/vasp5n/vasp.5.2/vasp
cp OUTCAR OUTCAR_soc110
cp OSZICAR OSZICAR_soc110

cp INCAR_soc111 INCAR
custodian_static.py
#mpirun -n $nCores /opt/vasp5n/vasp.5.2/vasp
cp OUTCAR OUTCAR_soc111
cp OSZICAR OSZICAR_soc111

clearFiles.sh
