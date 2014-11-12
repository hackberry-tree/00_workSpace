#!/bin/bash
#===========================================================
# explanation;
# make disordered structure for CVM from ATAT
#
# comment;
# need modification for adopting multi-sublattice CVM,
# and to be expected to marge dis-ord 
#
# 26/Mar/2011
# 31/Mar/2012 simplify
#===========================================================

home_directory=$PWD
str_num=$1
nph=$2

phase=`sed -n 3,3p input`
element=(`sed -n 1,1p input`)
set -- `sed -n 1,1p input`
num_ary=`expr $# - 1`

# neel1: number of atom positions
i=0
while [ ${i} -le ${num_ary} ]
do
   neel1[${i}]=0
   i=`expr ${i} + 1`
done

file_input=$str_num"/"str.out
file_output="CVM/"$phase"_dis.str"

#	crystal system
if [ $phase = "fcc" ]; then
           cph="FCC"
           brav="F"
elif [ $phase = "bcc" ]; then
           cph="BCC"
           brav="I"
elif [ $phase = "hcp" ]; then
           cph="HCP"
           brav="P"
else
           cph="PRM"
           brav="P"
fi

a=`sed -n 1,1p $file_input`
b=`sed -n 2,2p $file_input`
c=`sed -n 3,3p $file_input`

pname=`cut -d " " -f 1 pname.tmp`
echo "${cph}0*${pname} BRAV=${brav} DISORD=${cph}0 NPH="${nph} >> $file_output
echo "        A_CART="$a >> $file_output
echo "        B_CART="$b >> $file_output
echo "        C_CART="$c >> $file_output

#===========================================================
#  reading coordination of atoms
#	neel : total number of atom positions,
#	neel#: number of "#" atom positions
#	eel1 : name of element in the input list
#===========================================================
neel=0
nnf=7
eel1=`sed -n $nnf,$nnf\p $file_input`

_IFS="$IFS"
IFS=" "
set -- $eel1
pos_x=$1
pos_y=$2
pos_z=$3
name_elem=$4
IFS="$_IFS"

A=(A B C D E F G H I J K L M N O P Q R S T U V W X Y Z)
while [ $name_elem ]; do
     i=0
     while [ ${i} -le ${num_ary} ]
     do
        if [ "$name_elem" = "${element[${i}]}" ]; then
           echo "        ""ATOM"$(( neel + 1 ))"="${A[${i}]}" "$pos_x" "$pos_y" "$pos_z >> $file_output
           neel=`expr $neel + 1`
        fi
        i=`expr ${i} + 1`
     done

# reading next coordination of atoms
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


