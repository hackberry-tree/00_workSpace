#!/bin/bash
FILEIN=POSCAR.ini
FILEOU=CONTCAR

MIN=1
MAX=4 

rm SUMMARY.pos
for I in `seq $MIN $MAX`;
do
 echo ' ----- position ABCDE:' $I >>./SUMMARY.pos
 cd $I/
 NB=`cat -n $FILEIN | grep 'A01' | awk '{print $1}' `
 var=`head -$NB $FILEOU | tail -1`
 echo $var >>../SUMMARY.pos
 NB=`cat -n $FILEIN | grep 'B01' | awk '{print $1}' `
 var=`head -$NB $FILEOU | tail -1`
 echo $var >>../SUMMARY.pos
 cd ../
done
exit 0;
