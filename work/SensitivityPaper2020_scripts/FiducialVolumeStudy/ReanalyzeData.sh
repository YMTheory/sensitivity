#!/bin/bash

FIDVOL=$1
DATADIR=/p/lustre2/lenardo1/sensitivity_output/Apr10_SensitivityVsFiducialVolume_D024/

for filename in $(ls $DATADIR*$FIDVOL*[0-9].h5)
do
echo $filename
python ../Reanalyze90PercentLimitWithNewCriticalLambda.py $filename $(ls CriticalLambdaCurves/*$FIDVOL*) $DATADIR  
done 
