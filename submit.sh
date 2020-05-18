#!/bin/bash

if [ -z "$1" ]
then
    echo "Usage: ./`basename $0` comment njobs[=1] nevents[=10] bmin[=0] bmin[=5]"
    echo ""
    echo "bmin and bmax correspond the lower and upper limits of the centrality bins:"
    echo ""
    echo "    cent class | bmin | bmax  "
    echo "   ---------------------------"
    echo "     5 - 10 %  |  0   |  1    "
    echo "    10 - 20 %  |  1   |  2    "
    echo "    20 - 30 %  |  2   |  3    "
    echo "    30 - 40 %  |  3   |  4    "
    echo "    40 - 50 %  |  4   |  5    "
    echo ""
    exit 0
fi

if [ -z "$2" ]
then
    njobs=1
else
    njobs=$2
fi

if [ -z "$3" ]
then
    nevents=10
else
    nevents=$3
fi

if [ -z "$4" ]
then
    bmin=0
else
    bmin=$4
fi

if [ -z "$5" ]
then
    bmax=5
else
    bmax=$5
fi


for (( i=1; i<=$njobs; i++ ))
do
    outputdir=run_${1}_job$i
    mkdir $outputdir
    mkdir ${outputdir}/logs
    sbatch -o ${outputdir}/logs/log$i -e ${outputdir}/logs/errout$i -J pHI -n 1 run $i $nevents $outputdir $bmin $bmax
    sleep 1
done
