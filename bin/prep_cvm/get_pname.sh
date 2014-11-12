#!/bin/bash

#-----------------------------------------------------------
# explanation;
# analyse str.out as following input and print consentration
#
# usage; 
# get_pname.sh # (#: structure number at the directory)
#
# example;
# ./get_pname.sh 10
#            -> in $PWD/pname.tmp
#
# 30/Jan/2012
# 31/Mar/2012 simplify
# 13/Feb/2013 reconstract
#-----------------------------------------------------------

home_directory=$PWD

if [ -z $1 ]; then
   echo "FATAL"
   echo "can't find directory"
   echo "enter the directory after get_pname.sh"
   exit
else
   str_num=$1
fi

if [ ! -e "$home_directory/input" ]; then
   echo "FATAL"
   echo "can't find input in this directory"
   echo "make input before running get_pname.sh"
   exit
else
   element=(`sed -n 1,1\p input`)
   set -- `sed -n 1,1\p input`
   num_ary=`expr $# - 1`
fi

# files
file_input=$str_num"/"str.out
file_out="$home_directory/pname.tmp"

i=0
while [ ${i} -le ${num_ary} ]
do
# neel1: number of atom positions
   neel1[${i}]=0
   i=`expr ${i} + 1`
done

# nnf : start line reading str.out
nnf=7

# neel : total number of atom positions,
neel=0

# eel1 : name of element in the input list
eel1=`sed -n $nnf,$nnf\p $file_input`

_IFS="$IFS"
IFS=" "
set -- $eel1
pos_x=$1
pos_y=$2
pos_z=$3
name_elem=$4
IFS="$_IFS"

while [ $name_elem ]; do   
   i=0
   while [ ${i} -le ${num_ary} ]
   do
      if [ "$name_elem" = "${element[${i}]}" ]; then
         neel1[${i}]=`expr ${neel1[${i}]} + 1`
      fi
      i=`expr ${i} + 1`
   done


   if [ "$name_elem" != "Vac" ]; then
      neel=`expr $neel + 1`
   fi

#-----------------------------------------------------------
#       reading next coordination of atoms
#-----------------------------------------------------------

   nnf=`expr $nnf + 1`
   eel1=`sed -n $nnf,$nnf\p $file_input`
   _IFS="$IFS"
   IFS=" "
   set -- $eel1
   pos_x=$1
   pos_y=$2
   pos_z=$3
   name_elem=$4
done

#-----------------------------------------------------------
#       output pname
#-----------------------------------------------------------
i=0
pname=""
A=(A B C D E F G H I J K L M N O P Q R S T U V W X Y Z)
#echo ${num_ary}
while [ ${i} -le ${num_ary} ]
do
   pname=${pname}${A[${i}]}${neel1[${i}]}
   i=`expr ${i} + 1`
done

echo "$pname $neel" > $file_out
#

