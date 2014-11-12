#!/bin/bash
#===========================================================
# explanation;
# file preparation from ATAT to CVM
#
# usage; ./prep_cvm.sh
#
# call; 
# make_input.sh exte_vasp.sh exte_fit.sh
# constr_dis.sh constr_ord.sh get_pname.sh
# get_atom.sh  
#
# CVM template           -> in $PWD/CVM  directory
#
# 27/Apr/2011
#  4/Dec/2011   structure numbers are extracted from fit.out
# 30/Jan/2012   improve to deal with multicomponent
# 13/Feb/2013   reconstract
# 06/Nov/2014	
#===========================================================

home_directory=$PWD

if [ -d "$home_directory/CVM" ]; then
   echo Warning: Directory CVM found, replace it...
   rm -rf $home_directory/CVM
   mkdir CVM
else
   mkdir CVM
fi

if [ -e input ]; then
   echo Warning: File input found, replace it...
   rm input
fi

#===========================================================
#	reading lat.in
#===========================================================
echo ""
echo "[Information of ATAT]"

elements_all=`get_atom.sh -a`
set -- $elements_all
tot_num_elem=$#
echo All elements list:          $elements_all

elements_fixed=`get_atom.sh -f`
set -- $elements_fixed
num_elem=$#
if [ $num_elem != 0 ]; then
   echo Perfectly ordered elements list in ATAT, ${elements_fixed} 
fi
elements=`get_atom.sh -v`
set -- $elements
num_elem=$#
echo Elements considered in CVM: $elements

#===========================================================
#	creating input file
#===========================================================

make_input.sh

#===========================================================
#	preparation for energies.txt
#===========================================================

phase=`sed -n 3,3p input`
nph_base=`sed -n 4,4p input`
nph=`sed -n 5,5p input`
str_num_list=`sed -n 6,6p input`
base_str_num_list=`sed -n 7,7p input`
exte_input=`sed -n 8,8p input`

if [ v = $exte_input ];then
   ext_command="exte_vasp.sh"
elif [ f = $exte_input ];then
   ext_command="exte_fit.sh"
fi

echo ""
echo "[Start creating energies & structure file]"
echo "Base structure data retrieved from: "#/str.out
echo "Structure data retrieved from:      "#/str.out
echo "Energy data retrieved from:         "#/OSZICAR.static

for j in $str_num_list
do
   str=nobase
   for i in $base_str_num_list
   do
      if [ "${j}" = "${i}" ]; then
         str=base
      fi
   done

   if [ "$str" = "base" ]; then
      #===========================================================
      #		preparation for disordered structure
      #===========================================================
      echo "STR = ${j}, NPH = $nph_base"
      get_pname.sh ${j}
      constr_dis.sh ${j} ${nph_base}
      $ext_command ${j}
   else
      #===========================================================
      # 	preparation for ordered structure
      #===========================================================
      echo "STR = ${j}, NPH = $nph"
      get_pname.sh ${j}
      constr_ord.sh ${j} $nph
      $ext_command ${j}
      nph=`expr $nph + 1`
   fi

done

cat CVM/${phase}_dis.str CVM/${phase}_ord.str >CVM/${phase}.str

#==================================
#   fix energies.txt
#==================================
for i in $base_str_num_list; do
    pre=`grep "[A-Z]$i\*" CVM/energies.txt | sed "s/\*/\./g"`
    post=`grep "[A-Z]$i\*" CVM/energies.txt | sed "s/$i\*/0\*/g"`
    sed "s/$pre/$post/" CVM/energies.txt > CVM/energies.txt_tmp
    mv CVM/energies.txt_tmp CVM/energies.txt
done

#   fix num atom
sed "s/A/  /g" CVM/energies.txt | sed "s/B/  /g" | awk '{print $2+$3}' > atoms_tmp
paste CVM/energies.txt atoms_tmp | awk '{print $1, $2, $3, $5}' > energies_tmp
mv energies_tmp  CVM/energies.txt


#==================================
#   make IN.CVM
#==================================
strID=`grep "[A-Z]0\*" CVM/energies.txt | head -n 1 | sed "s/\s.*//g"`
pdir=`which prep_cvm.sh | sed 's/prep_cvm.sh//g'`
eTOT=`echo $strID | sed "s/[A-Z]*0\*//g" | sed "s/[0-9]\+\$/=\!/g" | sed "s/^/ETOT/g" | sed "s/[0-9]\+/=\!  ETOT/g"`
numSp=`echo $strID | sed "s/[A-Z]*0\*//g" | sed "s/[A-Z][0-9]\+/ A/g" | wc -w`

cMin=
cStp=
cSpn=
count=$numSp
while test $count -gt 2; do
 cMin="${cMin} 0.01"
 cStp="${cStp} 0"
 cSpn="${cSpn} 0.25"
 count=`expr $count - 1`
done

cat $pdir/IN.CVM | sed "s/__StrID__/$strID/g" | sed "s/__ETOT__/$eTOT/g" | sed "s/__STR__/$phase/g" | sed "s/__NumSp__/$numSp/g" | sed "s/__C__/$cMin/g" | sed "s/__CSTP__/$cStp/g" | sed "s/__CSPN__/$cSpn/g" > CVM/IN.CVM

cp $pdir/def.txt CVM/
