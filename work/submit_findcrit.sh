!/bin/bash

# Author: Samuele Sangiorgio <samuele@llnl.gov>


## SETTINGS

SIGNALLIST1=$(seq 0.0 3 15.0)
SIGNALLIST2=$(seq 7.0 0.3 10.) 


# number of runs (seeds) for each signal
SEEDS=$(seq 1 100)

# number of toys
NTOYS=1000

# Max number of jobs allowed in the queue
MAXJOBS=200

# LC username
LCUSER=stiegler

# Where all log files are saved (a folder for each run number is created). Default is current directory
LOGSDIR=$(pwd)/../logs

# Walltime
WALLTIME=6:00:00

# BATCH QUEUE
QUEUE=pbatch

#MAGICNUMDIR=/g/g19/brodsky3/nexo/scratch/sens_recalc/results/zerobkg/
OUTPUTDIR=/p/lscratche/nexouser/sens_recalc/results/

ROOTTABLE=/p/lscratche/nexouser/sens_recalc/tables/Summary_test_2.root

SCRIPTDIR=/p/lscratche/nexouser/sens_recalc/work

for SEED in $SEEDS
do
#    for SIGNAL in "${SIGNALLIST1[@]}" $SIGNALLIST2
#   for SIGNAL in $SIGNALLIST1 $SIGNALLIST2
#    do
	SIGNAL=0
        # Talk to the user
        echo ""
        echo "-----------------------------------------------------------------"
        echo "Starting on seed number: $SEED and signal $SIGNAL"

#       # check if run has been already submitted
        if grep -q -x "${SIGNAL}-${SEED}" ${LOGSDIR}/submitted.txt ; then
            echo "ERROR: Run number already in the queue or done"
            continue
        fi
        echo 
        DIR=${LOGSDIR}/$SIGNAL-$SEED
        echo "Logs Output Directory: ${DIR}"
        [ ! -d ${DIR} ] && mkdir $DIR # Create folder if doesn't exist

        # now enter the run folder
        cd ${DIR}

        # clean up the logs 
        rm slurm-* *.core run.* runtime-*  &>/dev/null

        # create the run script
        echo "#!/bin/bash" > run.sh
        echo "echo START TIME: [\`date\`]" >> run.sh
        echo "START=\`date +%s\`"  >> run.sh
        echo "cd ${SCRIPTDIR}" >> run.sh
        # THIS IS FOR THE FIRST STEP TO CALCULATE LAMBDA CRITICAL
        echo "python RunSensitivity.py -n ${NTOYS} -r 1 -b -y 10.0 -c ${SIGNAL} -s ${SEED} -d ${OUTPUTDIR} -o allbkgs -t ${ROOTTABLE} --ssfrac-improvement 1.0 --rn222-rate-correction 1.00 --turn-off-groups GhostComponents"  >> run.sh
	#echo "python RunSensitivity.py -n ${NTOYS} -r 1 -b -y 10.0 -c ${SIGNAL} -s ${SEED} -d ${OUTPUTDIR} -o allbkgs -t ${ROOTTABLE} --ssfrac-improvement 1.0 --rn222-rate-correction 1.00 --turn-off-groups GhostComponents"  >> run.sh
        # THIS IS FOR THE SECOND STEP OF EXTRACTING UPPER CASES
        #echo "python MakeMagicNumberTable.py -n ${NTOYS} -r 0 -y 10.0 -c ${SIGNAL} -s ${SEED} -d ${OUTPUTDIR} -o recalc_test -t ${ROOTTABLE} -m 0 --ssfrac-improvement 1.0 --rn222-rate-correction 1.00 -M ${MAGICNUMDIR}"  >> run.sh
        echo "cd -" >> run.sh
        echo "END=\`date +%s\`"  >> run.sh
        echo "RUNTIME=\$((END-START))"  >> run.sh
        echo "echo RUNTIME: \$RUNTIME seconds | tee runtime-${SIGNAL}-${SEED}.txt"   >> run.sh
        echo "echo STOP TIME: [\`date\`]" >> run.sh
        chmod 700 run.sh

        # if more than MAXJOBS jobs in queue, wait before submitting
        while [ $(squeue -u ${LCUSER} | wc -l) -gt $MAXJOBS ]
        do
            echo "Waiting for room in the queue"        
            sleep 60    # in seconds
        done

        # execute or send to the queue
        #./run.sh 2>run.err | tee run.out
        msub -A afqn -N ${SIGNAL}-${SEED} -l walltime=${WALLTIME} -q ${QUEUE} -V -o slurm-${SIGNAL}-${SEED}.out -e slurm-${SIGNAL}-${SEED}.err run.sh
                
        echo ${SIGNAL}-${SEED} >> ${LOGSDIR}/submitted.txt

    #done
done
