#!/bin/bash

cp KPOINTS_relax_reduced KPOINTS
cp INCAR_volumeE INCAR
custodian_relax.py
#mpirun -n $nCores /opt/vasp5/vasp.5.2/vasp
#judgeRX.py

cp INCAR_cell INCAR
custodian_relax.py
#mpirun -n $nCores /opt/vasp5/vasp.5.2/vasp
#judgeRX.py

clearFiles.sh
