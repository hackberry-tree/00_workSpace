#!/bin/bash 
i=$(python /home/enoki/Dropbox/00_scripts/bin/batch_run.py $1)
while $i
do
     i=$(python /home/enoki/Dropbox/00_scripts/bin/batch_run.py $1)
     sleep 300
done
