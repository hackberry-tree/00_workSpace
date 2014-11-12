#!/bin/bash
cp KPOINTS_relax KPOINTS
cp INCAR_relax INCAR

custodian_relax.py

clearFiles.sh
