#!/bin/bash

if [ "$1" == "-h" ]
then
    echo "Usage: `basename $0` comment njobs[=1]"
    exit 0
fi

if [ -z "$1" ]
then
    echo "Please give a comment to make this run unique (check `basename $0` -h for help)"
    exit 0
fi

if [ -z "$2" ]
then
    njobs=1
else
    njobs=$2
fi

outputdir=run_${1}
mkdir $outputdir
mkdir ${outputdir}/logs

for (( i=1; i<=$njobs; i++ ))
do
    sbatch -o ${outputdir}/logs/log$i -e ${outputdir}/logs/errout$i -J pHI -n 1 run
    sleep 1
done
