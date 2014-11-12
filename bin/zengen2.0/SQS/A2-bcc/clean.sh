#!/bin/bash
#

for DIR in *; do
  if [ -d $DIR ]; then
   echo "directory = $DIR"
   cd $DIR
   rm WAVECAR
   rm CHG*
   rm PCDAT PROCAR vasprun.xml XDATCAR  EIG*
   rm OUTCAR.*
   cd ..
  fi
done
exit 0;
