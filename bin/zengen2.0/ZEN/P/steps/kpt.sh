#!/bin/bash

for a in `seq 2 6` ;
do

# CHECK THE RATIO !
b=$(echo "scale=2; $a*3.58916/1" | bc | awk '{printf("%d\n",$0+=$0<0?-0.5:0.5)}')
c=$(echo "scale=2; $a*1.90647/1" | bc | awk '{printf("%d\n",$0+=$0<0?-0.5:0.5)}')

echo "****************************** ka= $a"
echo "********* kb= $b"
echo "********* kc= $c"

mkdir $a
cd $a
cp ../POTCAR .
cp ../POSCAR .
cp ../INCAR .
cp ../WAV* .
cp ../CHG* .

cat >KPOINTS <<!
k mesh
 0
Gamma           
$a $b $c   
 0  0  0
!

# RUN VASP
llsubmit ../job

cd ..

done

exit 0;

