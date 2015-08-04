#!/bin/bash
i=$(python /home/enoki/Dropbox/Codes/00_workSpace/bin/batch_run.py $1)
while $i
do
     i=$(python /home/enoki/Dropbox/Codes/00_workSpace/bin/batch_run.py $1)
     sleep 300
done
