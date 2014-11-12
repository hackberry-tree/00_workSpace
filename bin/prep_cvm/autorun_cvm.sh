#!/bin/bash
icvm
Diso=`grep "[A-Z]0\*" energies.txt | head -1 | sed "s/\*.*//"`
while test `grep "Same" log.txt | grep -v $Diso | wc -l` -ne 0; do
    python ../debug_cvm.py
    rm cwm.txt cmt.txt
    icvm
done
