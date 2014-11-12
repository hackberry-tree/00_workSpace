#!/bin/bash
VER=5 # vasp 4 or 5
echo -n "Because of CONTCAR files format: Which Version of Vasp used (4 or 5) ?"
read -e VER

FILEIN=POSCAR.ini
FILEOU=CONTCAR
echo -n "From MIN to MAX Compounds ?"
read -e MIN
read -e MAX

rm SUMMARY.pos
echo "script for internal positions of chi-phase - VASP " $VER 
VER=`expr $VER - 4`

for I in `seq $MIN $MAX`;
do
 echo ' ----- position ABCDE:' $I >>./SUMMARY.pos
 cd $I/
 NB=`cat -n $FILEIN | grep 'A01' | awk '{print $1}' `
 NB=`expr $NB + $VER`
 var=`head -$NB $FILEOU | tail -1`
 echo $var >>../SUMMARY.pos
 NB=`cat -n $FILEIN | grep 'B01' | awk '{print $1}' `
 NB=`expr $NB + $VER`
 var=`head -$NB $FILEOU | tail -1`
 echo $var >>../SUMMARY.pos
 NB=`cat -n $FILEIN | grep 'C01' | awk '{print $1}' `
 NB=`expr $NB + $VER`
 var=`head -$NB $FILEOU | tail -1`
 echo $var >>../SUMMARY.pos
 NB=`cat -n $FILEIN | grep 'D01' | awk '{print $1}' `
 NB=`expr $NB + $VER`
 var=`head -$NB $FILEOU | tail -1`
 echo $var >>../SUMMARY.pos
 NB=`cat -n $FILEIN | grep 'E01' | awk '{print $1}' `
 NB=`expr $NB + $VER`
 var=`head -$NB $FILEOU | tail -1`
 echo $var >>../SUMMARY.pos
 NB=`cat -n $FILEIN | grep 'F01' | awk '{print $1}' `
 NB=`expr $NB + $VER`
 var=`head -$NB $FILEOU | tail -1`
 echo $var >>../SUMMARY.pos
 NB=`cat -n $FILEIN | grep 'G01' | awk '{print $1}' `
 NB=`expr $NB + $VER`
 var=`head -$NB $FILEOU | tail -1`
 echo $var >>../SUMMARY.pos
 NB=`cat -n $FILEIN | grep 'H01' | awk '{print $1}' `
 NB=`expr $NB + $VER`
 var=`head -$NB $FILEOU | tail -1`
 echo $var >>../SUMMARY.pos
 NB=`cat -n $FILEIN | grep 'I01' | awk '{print $1}' `
 NB=`expr $NB + $VER`
 var=`head -$NB $FILEOU | tail -1`
 echo $var >>../SUMMARY.pos
 NB=`cat -n $FILEIN | grep 'J01' | awk '{print $1}' `
 NB=`expr $NB + $VER`
 var=`head -$NB $FILEOU | tail -1`
 echo $var >>../SUMMARY.pos
 NB=`cat -n $FILEIN | grep 'K01' | awk '{print $1}' `
 NB=`expr $NB + $VER`
 var=`head -$NB $FILEOU | tail -1`
 echo $var >>../SUMMARY.pos
 NB=`cat -n $FILEIN | grep 'L01' | awk '{print $1}' `
 NB=`expr $NB + $VER`
 var=`head -$NB $FILEOU | tail -1`
 echo $var >>../SUMMARY.pos
 cd ../
done
exit 0;
