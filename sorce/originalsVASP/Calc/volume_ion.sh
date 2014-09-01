#!/bin/bash

#optimize volume
cp KPOINTS_relax_reduced KPOINTS
cp INCAR_volume INCAR
custodian_relax.py

mkdir vol_errors
mv error.* vol_errors

#optimize ion
cp INCAR_ion INCAR
custodian_relax1.py

clearFiles.sh
