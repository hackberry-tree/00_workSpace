#!/bin/bash
cp KPOINTS_soc KPOINTS
cp INCAR_dos INCAR

custodian_static.py
#mpirun -n $nCores /opt/vasp5/vasp.5.2/vasp
cp DOSCAR DOSCAR_polarized

cp KPOINTS_bandTet KPOINTS ###保留
cp INCAR_band INCAR

custodian_static.py
#mpirun -n $nCores /opt/vasp5/vasp.5.2/vasp

cp OUTCAR OUTCAR_band
cp PROCAR PROCAR_band

