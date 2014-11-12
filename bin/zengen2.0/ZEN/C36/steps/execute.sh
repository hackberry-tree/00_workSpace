#!/bin/bash
# ver3 - june 2011

 ################
 ### Parameters #
 ################
NbCPU=4 # Number of CPU allowed
# Which compound
CMIN=1  # Minimal Compound number 
CMAX=2  # Maximal Compound number 
# Which STEP ?
# 0:         initialization
# 1,2,3,4... 
SMIN=0  # Minimal relaxation step 
SMAX=0  # Maximal relaxation step 
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
 cp POSCAR CONTCAR
 else
   cp CONTCAR POSCAR
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

