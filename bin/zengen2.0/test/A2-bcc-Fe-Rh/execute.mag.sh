#!/bin/bash
# spin-polarized sqs     

 ################
 ### Parameters #
 ################

 # Number of cores ?
NbCPU=12 # Number of CPU allowed

 # Which STEP ?
 # With spin polarized calculation : 2 steps available
SMIN=1  # Minimal relaxation step 
SMAX=2  # Maximal relaxation step  
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
 if [ $STEP -eq 1 ] ; then
 # new POSCAR
 cp CONTCAR.1 POSCAR
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
 cp CONTCAR CONTCAR.M$STEP
 cp OUTCAR OUTCAR.M$STEP
 cp OSZICAR OSZICAR.M$STEP
 E=`tail -1 OSZICAR`
 echo $DIR $E >>../SUMMARY.M$STEP
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
echo "Hello, "$USER", spin polarized calculation." 
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
 echo "relaxstep = M$STEP *******************************"
 echo "directory = $DIR *********************************"
 cd $DIR
 preparposcar 
 # initialisation
 cp ../INCAR.M$STEP INCAR
 cp KPOINTS.$STEP KPOINTS
 #
 execution
 savefiles
 cd ..
##############################
  fi
done
 # summary            
 if [ $STEP -eq 1 ] ; then
 perl sum.pl
 mv sum.out sum-mag.vol.out        
 fi
 if [ $STEP -eq 2 ] ; then
 perl sum.pl
 mv sum.out sum-mag.full.out        
 fi
done
##############################
echo "*********************************"
exit 0;
###################################
###################################

