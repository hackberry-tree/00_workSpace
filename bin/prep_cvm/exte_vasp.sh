#!/bin/bash
#-----------------------------------------------------------
# retrive VASP energies
# usage; exte_vasp.sh #
#                              (#: structure number at the directory)
#
# results of VASP        -> in $PWD/CVM  directory
#            cat vasp_h_cmat.txt vasp_cmat.txt> cmat.txt
#            cat vasp_h_energy.txt vasp_energy.txt > energy.txt
#           (i.e., for example, cat ) 
#
# 24/Mar/2011
# 12/Feb/2013 simplify
#-----------------------------------------------------------

home_directory=$PWD
str_num=$1

element=(`get_atom.sh -a`)
set -- `get_atom.sh -a`
num_ary=`expr $# - 1`

# neel1: number of atom positions
i=0
while [ ${i} -le ${num_ary} ]
do
   neel1[${i}]=0
   i=`expr ${i} + 1`
done

# files
filename_input=$home_directory"/"$str_num"/"OSZICAR.static
file_input=$str_num"/"str.out
file_vasp_out=$home_directory/CVM/energies.txt

phase=`sed -n 3,3p input`
if [ $phase = "fcc" ]; then
           cph="FCC"
           ccph="FCC"
elif [ $phase = "bcc" ]; then
           cph="BCC"
           ccph="BCC"
elif [ $phase = "hcp" ]; then
           cph="HCP"
           ccph="HCP"
else
           cph="PRM"
           ccph="PRM"
fi

#-----------------------------------------------------------
# neel : total number of atom positions,
# neel1: number of "1" atom positions
# neel2: number of "2" atom positions
# eel1 : name of element in the input list
#-----------------------------------------------------------

nnf=7
neel=0
ref_energy=`sed -n 2,2p input`
eel1=`sed -n $nnf,$nnf\p $file_input`

_IFS="$IFS"
IFS=" "
set -- $eel1
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
     name_elem=$4
done

#-----------------------------------------------------------
#       output
#-----------------------------------------------------------
pname=`cut -d " " -f 1 pname.tmp`

i=0
ref_sum=0
for ref_ene1 in $ref_energy
do
#   echo $ref_ene1
   ref_sum=`echo "scale=10;  $ref_sum + $ref_ene1 * ${neel1[${i}]}" | bc`
   i=`expr ${i} + 1`
done

#echo $ref_sum
if [ $ref_sum = 0 ]; then
   energy_vasp=`grep F */OSZICAR.static | grep ^${str_num}/ |cut -d " " -f 6`
else
   energy_vasp=`grep F */OSZICAR.static | grep ^${str_num}/ |cut -d " " -f 6| sed s/E+0/*10^/g`
   energy_vasp=`echo "scale=10;  $energy_vasp - $ref_sum" | bc`
fi

#energy_out=$cph$str_num"*"$pname" 1 "$energy_vasp" "$neel
energy_out=$cph$str_num"*"$pname" 1 "$energy_vasp" "$neel
echo $energy_out >> $file_vasp_out

IFS="$_IFS"


 


