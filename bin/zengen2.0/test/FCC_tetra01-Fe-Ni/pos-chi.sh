#!/bin/bash
# version 1.3 - jc.crivello - jan 2012

VER=5 # vasp 4 or 5
echo -n "Because of CONTCAR files format: Which Version of Vasp used (4 or 5) ?"
read -e VER

FILEIN=POSCAR.ini
FILEOU=CONTCAR
echo -n "From MIN to MAX Compounds ?"
read -e MIN
read -e MAX


##########################################################
rm pos.out
echo "script for internal positions of chi-phase - VASP " $VER 
VER=`expr $VER - 4`
for I in `seq $MIN $MAX`;
do
 echo ' ----- position ABCD in folder:' $I 
 echo ' ----- position ABCD:' $I >>./pos.out
 cd $I/
 NB=`cat -n $FILEIN | grep 'A01' | awk '{print $1}' `
 NB=`expr $NB + $VER`
 var=`head -$NB $FILEOU | tail -1`
 echo $var >>../pos.out
 NB=`cat -n $FILEIN | grep 'B01' | awk '{print $1}' `
 NB=`expr $NB + $VER`
 var=`head -$NB $FILEOU | tail -1`
 echo $var >>../pos.out
 NB=`cat -n $FILEIN | grep 'C01' | awk '{print $1}' `
 NB=`expr $NB + $VER`
 var=`head -$NB $FILEOU | tail -1`
 echo $var >>../pos.out
 NB=`cat -n $FILEIN | grep 'D01' | awk '{print $1}' `
 NB=`expr $NB + $VER`
 var=`head -$NB $FILEOU | tail -1`
 echo $var >>../pos.out
 cd ../
done
exit 0;
