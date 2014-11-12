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
file_output="CVM/"$phase"_ord.str"

#	crystal system
if [ $phase = "fcc" ]; then
           cph="FCC"
elif [ $phase = "bcc" ]; then
           cph="BCC"
elif [ $phase = "hcp" ]; then
           cph="HCP"
else
           cph="PRM"
fi

# conversion from abc & uvw combination to a,b,c_cart
           a=`sed -n 1,1p $file_input`
           b=`sed -n 2,2p $file_input`
           c=`sed -n 3,3p $file_input`
           u=`sed -n 4,4p $file_input`
           v=`sed -n 5,5p $file_input`
           w=`sed -n 6,6p $file_input`

           _IFS="$IFS"
           IFS=" "
           set -- $a
           ax=$1
           ay=$2
           az=$3
           set -- $b
           bx=$1
           by=$2
           bz=$3
           set -- $c
           cx=$1
           cy=$2
           cz=$3
           set -- $u
           ua=$1
           ub=$2
           uc=$3
           set -- $v
           va=$1
           vb=$2
           vc=$3
           set -- $w
           wa=$1
           wb=$2
           wc=$3
           IFS="$_IFS"

           a_cart_1=`echo "scale=10;  $ua*$ax+$ub*$bx+$uc*$cx" | bc`
           a_cart_2=`echo "scale=10;  $ua*$ay+$ub*$by+$uc*$cy" | bc`
           a_cart_3=`echo "scale=10;  $ua*$az+$ub*$bz+$uc*$cz" | bc`
           a_cart=$a_cart_1" "$a_cart_2" "$a_cart_3
           
           b_cart_1=`echo "scale=10;  $va*$ax+$vb*$bx+$vc*$cx" | bc`
           b_cart_2=`echo "scale=10;  $va*$ay+$vb*$by+$vc*$cy" | bc`
           b_cart_3=`echo "scale=10;  $va*$az+$vb*$bz+$vc*$cz" | bc`
           b_cart=$b_cart_1" "$b_cart_2" "$b_cart_3
           
           c_cart_1=`echo "scale=10;  $wa*$ax+$wb*$bx+$wc*$cx" | bc`
           c_cart_2=`echo "scale=10;  $wa*$ay+$wb*$by+$wc*$cy" | bc`
           c_cart_3=`echo "scale=10;  $wa*$az+$wb*$bz+$wc*$cz" | bc`
           c_cart=$c_cart_1" "$c_cart_2" "$c_cart_3

pname=`cut -d " " -f 1 pname.tmp`
echo "${cph}${str_num}*${pname} BRAV=P DISORD=${cph}0 NPH="${nph} >> $file_output
echo "        A_CART="$a_cart >> $file_output
echo "        B_CART="$b_cart >> $file_output
echo "        C_CART="$c_cart >> $file_output

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

#       reading next coordination of atoms
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
