#!/bin/bash
# Script to add the POTCAR file of an element x.
# ver_1.9 - may 2013

potdir=/home/enoki/vasp/vasp_potential/vasp/potpaw_PBE
#ex: potdir=/home/jcc/POT/potpaw_PBE

if [ $# -lt 1 ]; then
        echo "Usage: potcar [Element-name]"
        echo "ex: potcar Mn_pv ..."
        exit 1
fi

element=$1

#gzip -d $potdir/$element/POTCAR.Z
cp $potdir/$element/POTCAR POTCAR.$element

exit 0

