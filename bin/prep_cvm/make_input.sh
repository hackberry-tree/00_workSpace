#!/bin/bash
#===========================================================
#
# file preparation from ATAT to CVM
#
# usage; ./make_input.sh
#
# in $PWD/input
#
# 12/Feb/2013   sepalated ed by prep_cvm.sh
#===========================================================

home_directory=$PWD

echo ""
echo "[Create input file for the preparation]"

elements=""
#echo ""
echo "input elements assigned for A B C... (delimiter is space) "
read elements

echo ""
echo "input reference energy for $elements_all (delimiter is space) "
read ref_energy

echo ""
echo "input base structure (fcc, bcc, hcp and others)"
read phase

echo ""
echo "input NPH for base structure (1-9)"
read nph_base

echo ""
echo "input NPH for first ordered structure (NPH>=10)"
read nph1

elements1=`get_atom.sh -v`
set -- $elements1
num_elem=$#

if [ $num_elem = 2 ]; then
   echo ""
   # echo "input structure numbers from fit.out "
   # str_num_list=`cut -d " " -f 6 fit.out`
   echo "input structure numbers from exist */energy "
   str_num_list=`find */energy | sed "s/\/energy//g"`
   echo $str_num_list
else
   column_list=`expr $num_elem + 5`
   echo ""
   echo "input structure numbers from fit.out "
   str_num_list=`cut -d " " -f $column_list fit.out`
   echo $str_num_list
fi

echo "input a (accept) or structure numbers"
read man_str_num
if [ "$man_str_num" != "a" ]; then
   str_num_list=$man_str_num
fi

echo ""
echo "input base structure numbers (delimiter is space) "
read base_str_num_list

#echo ""
#echo "E_vasp(OSZICAR) or F_form(fit.out) (-> type v/f)"
#echo "Notice; fit.out is not available"
#read exte_input
exte_input="v"

echo $elements>>input
echo $ref_energy>>input
echo $phase >>input
echo $nph_base >>input
echo $nph1 >>input
echo $str_num_list >>input
echo $base_str_num_list >>input
echo $exte_input >>input

