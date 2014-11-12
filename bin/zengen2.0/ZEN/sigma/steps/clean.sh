#!/bin/bash
#
CMIN=3  # Minimal Compound number 
CMAX=32  # Maximal Compound number 

for DIR in `seq $CMIN $CMAX`;
 do
  echo "directory = $DIR"
#  cd $DIR
     rm $DIR/WAVECAR 
     rm $DIR/CHG*   
     rm $DIR/PCDAT $DIR/PROCAR $DIR/XDATCAR 
#  cd ..
done
exit 0;
~                  
