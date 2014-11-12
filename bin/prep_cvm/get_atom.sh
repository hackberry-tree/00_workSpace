#!/bin/bash
#-----------------------------------------------------------
# explanation;
# called by some scripts of structure conversion from ATAT to CVM
# read lat.in
# > display all species including Vacancy in the first line
# > second line shows F (Fixed) or V (Variable) consentration for each species
#
# comment;
# needs modification for adopting multi-sublattice CVM
# 
# 30/Jan/2012
# 20/Feb/2013
#-----------------------------------------------------------

home_directory=$PWD

while [ $# != 0 ]; do
  case $1
  in
    -a)
        showall=1;;
    -f)
        showfixed=1;;
    -v)
        showvariable=1;;
    -l)
        showlogical=1;;
    *)
    break;;
  esac
  shift
done

IFS=$'\n'
ln_file=(`cat lat.in`)

_IFS="$IFS"
IFS=" "

atom_ARRAY=("")
VF_ARRAY=("")
ln=0

for line in "${ln_file[@]}"; do
    ln=$(( ln + 1 ))

    set -- $line
    num_column=$#
    atoms=$4

    IFS=$','
    set -- $atoms
    num_elem=$#

    if [ $num_elem -gt 1 ]; then
        atom_var=V
    else
        atom_var=F
    fi


    _IFS="$IFS"
    IFS=" "

    while [ "$1" != "" ] ; do
        atom_list=not_exist 
# check first this atom is in the list or not 
        for (( i = 1; i < ${#atom_ARRAY[@]}; ++i )) ; do
            if [ $1 = ${atom_ARRAY[$i]}  ]; then
                atom_list=exist
            fi
        done
        if [ $atom_list = "not_exist" ]; then
            atom_ARRAY=("${atom_ARRAY[@]}" "$1")
            VF_ARRAY=("${VF_ARRAY[@]}" "$atom_var")
        fi
        shift 1
    done
done

if [ ! -z $showall ]; then
    echo ${atom_ARRAY[@]}
fi

if [ ! -z $showfixed ]; then
    fixedlist=""
    for (( i = 1; i < ${#atom_ARRAY[@]}; ++i )) ; do
        if [ ${VF_ARRAY[${i}]} = F ]; then
           fixedlist="$fixedlist ${atom_ARRAY[${i}]}" 
        fi
    done
    echo $fixedlist
fi

if [ ! -z $showvariable ]; then
    variablelist=""
    for (( i = 1; i < ${#atom_ARRAY[@]}; ++i )) ; do
        if [ ${VF_ARRAY[${i}]} = V ]; then
           variablelist="$variablelist ${atom_ARRAY[${i}]}" 
        fi
    done
    echo $variablelist
fi

if [ ! -z $showlogical ]; then
    echo ${VF_ARRAY[@]}
fi

