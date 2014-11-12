#!/bin/bash
# version 1.4 - jc.crivello - mar 2012

VER=4 # vasp 4 or 5
echo -n "Because of CONTCAR files format: Which Version of Vasp used (4 or 5) ?"
read -e VER

FILEIN=POSCAR.ini
FILEOU=CONTCAR
echo -n "From MIN to MAX Compounds ?"
read -e MIN
read -e MAX

##########################################################
rm pos.out
echo "script for internal positions of sigma-phase - VASP " $VER 
VER=`expr $VER - 4`
#for I in 02 03 04 05 06 07 08 09 `seq 10 $MAX`;
for I in `seq $MIN $MAX`;
do
 echo ' ----- position ABCDE in folder:' $I 
 echo ' ----- position ABCDE:' $I >>./pos.out
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
 NB=`cat -n $FILEIN | grep 'E01' | awk '{print $1}' `
 NB=`expr $NB + $VER`
 var=`head -$NB $FILEOU | tail -1`
 echo $var >>../pos.out
 cd ../
done
exit 0;
