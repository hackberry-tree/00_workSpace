#!/bin/bash

sh Calc/relax2stp.sh
vasp_voldep.py --5to10

dirs=`ls voldep/`
for dir in $dirs;
do
    cd voldep/$dir
    sh Calc/cell.sh
    sh Calc/static.sh
    cd ../../
done