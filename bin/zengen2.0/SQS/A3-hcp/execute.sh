#!/bin/bash
# ver4 - sqs january 2013 - 16 hcp 
# Structure from D. Shin          
# PRB 74 (2006) 024204      

 ################
 ### Parameters #
 ################

 # Number of cores ?
NbCPU=12 # Number of CPU allowed

 # Which STEP ?
SMIN=0  # Minimal relaxation step 
SMAX=6  # Maximal relaxation step  
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
 cp ../POTCAR .
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
 mpiexec -n $NbCPU ~/bin/vasp4m
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
##############################
for STEP in `seq $SMIN $SMAX`; do
for DIR in *; do
  if [ -d $DIR ]; then
##############################
 echo "relaxstep = $STEP *********************************"
 echo "directory = $DIR *********************************"
 cd $DIR
 preparposcar 
 # initialisation
 cp ../INCAR.$STEP INCAR
 cp KPOINTS.$STEP KPOINTS
 #
 execution
 savefiles
 cd ..
##############################
  fi
done
 # summary            
 if [ $STEP -eq 2 ] ; then
 perl sum.pl
 mv sum.out sum.vol.out        
 fi
 if [ $STEP -eq 5 ] ; then
 perl sum.pl
 mv sum.out sum.full.out        
 fi
done
##############################
echo "*********************************"
exit 0;
###################################
###################################

