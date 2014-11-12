#!/bin/bash
# ver4 - 2011

 ################
 ### Parameters #
 ################
NbCPU=4 # Number of CPU allowed
# Which compound
CMIN=1  # Minimal Compound number 
CMAX=3  # Maximal Compound number 
# Which STEP ?
# 0:         initialization
# 1,2,3,4:   draft calculations
# 5:         finalisation
SMIN=0  # Minimal relaxation step 
SMAX=5  # Maximal relaxation step 
# if single relaxation step
#SMAX=$SMIN 

###################################
 ################
 ### Functions ##
 ################
 # Prepa POSCAR
 ################
 preparposcar() 
 {
 if [ $STEP -eq 0 ] ; then
 # new POSCAR
 cp POSCAR.ini POSCAR
 else
 # diago CONTCAR
   CHA=$IFS
   cp CONTCAR TEST
   head -n 5 TEST > toto
   head -n 2 toto > toto1
   tail -n +3 toto > toto2
   tail -n +6 TEST > toto3
   IFS=$' '
   ligne1=($(sed -ne '1p' toto2))
   ligne2=($(sed -ne '2p' toto2))
   ligne3=($(sed -ne '3p' toto2))
   echo "Before:"
   echo ${ligne1[0]} ${ligne1[1]} ${ligne1[2]}
   echo ${ligne2[0]} ${ligne2[1]} ${ligne2[2]}
   echo ${ligne3[0]} ${ligne3[1]} ${ligne3[2]}
   AP=$(echo "scale=16; (${ligne1[1]})*1.0/1" | bc )
   AM=$(echo "scale=16; (${ligne1[1]})*-1.0/1" | bc )
   ligne1[0]=$AM
   ligne1[1]=$AP
   ligne1[2]=$AP 
   ligne2[0]=$AP 
   ligne2[1]=$AM 
   ligne2[2]=$AP 
   ligne3[0]=$AP 
   ligne3[1]=$AP 
   ligne3[2]=$AM 
   echo "After:"
   echo ${ligne1[0]} ${ligne1[1]} ${ligne1[2]}
   echo ${ligne2[0]} ${ligne2[1]} ${ligne2[2]}
   echo ${ligne3[0]} ${ligne3[1]} ${ligne3[2]}
   echo ${ligne1[0]} ${ligne1[1]} ${ligne1[2]} >toto2
   echo ${ligne2[0]} ${ligne2[1]} ${ligne2[2]} >>toto2
   echo ${ligne3[0]} ${ligne3[1]} ${ligne3[2]} >>toto2
   cat toto1 toto2 toto3 > POSCAR
   rm TEST toto toto1 toto2 toto3
   IFS=$CHA
 fi
 }
 ################
 # VASP execution 
 ################
 execution()
 {
 if [ $NbCPU -eq 1 ] ; then
 vasp4s
 else
 echo $PWD
 mpiexec -n $NbCPU vasp4m
 fi
 }
 ################
 # save files     
 ################
 savefiles()
 {
 cp CONTCAR CONTCAR.$STEP
 cp OUTCAR OUTCAR.$STEP
 cp OSZICAR OSZICAR.$STEP
 E=`tail -1 OSZICAR`
 echo $DIR $E >>../SUMMARY.$STEP
 T=`grep 'Elapsed time' 'OUTCAR'`
 echo $STEP $T >>TIME_SUMMARY.$DIR
 }
###################################


###################################
 #######
 # MAIN     
 #######
clear
ls
echo "*********************************"
echo "***********"
echo "Hello, "$USER"." 
echo "***********"
echo "ready to run on "$NbCPU" Nb of core(s) ?..."
echo "***********"
echo "From step #"$SMIN" to step #"$SMAX", ok?"
echo "***********"
echo "From compound #"$CMIN" to compound #"$CMAX", ok?"
echo "***********"
##############################
for STEP in `seq $SMIN $SMAX`; do
for DIR  in `seq $CMIN $CMAX`; do
##############################
 echo "relaxstep = $STEP *********************************"
 echo "directory = $DIR *********************************"
 cd $DIR
 preparposcar 
 # initialisation
 cp ../steps/INCAR.$STEP INCAR
 cp ../steps/KPOINTS.$STEP KPOINTS
 #
 execution
 savefiles
 cd ..
##############################
done
done
##############################
echo "*********************************"
exit 0;
###################################
###################################

