#!/bin/bash

#Make IBZKP
cp KPOINTS_soc KPOINTS
cp INCAR_ibzkp INCAR

mpirun -n $nCores /opt/vasp5n/vasp.5.2/vasp
cp OUTCAR OUTCAR_ibzkp
cp IBZKPT KPOINTS_ibzkp

vaspy -soc

#Polarized
cp INCAR_presoc INCAR
cp KPOINTS_ibzkp KPOINTS
mpirun -n $nCores /opt/vasp5/vasp.5.2/vasp
cp OSZICAR OSZICAR_presoc
cp OUTCAR OUTCAR_presoc

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
