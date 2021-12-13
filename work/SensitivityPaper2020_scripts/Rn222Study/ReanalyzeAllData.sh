#!/bin/bash

FILEDIR="/p/lustre2/lenardo1/sensitivity_output/Mar1_Rn222Study_merged-v11_OptimizedV1Binning_D024/"

for SCALE in 00.1x 00.3x 01.0x 03.0x 10.0x 100.0x
do

   for FILENAME in $(ls "$FILEDIR"*"$SCALE"*"h5")
   do 
      CRITICAL_LAMBDA_FILE=$(ls "CriticalLambdaCurves/critical_lambda_merged-v11_D024_"*"$SCALE"*"txt")
      python ../Reanalyze90PercentLimitWithNewCriticalLambda.py $FILENAME $CRITICAL_LAMBDA_FILE $FILEDIR
   done
        
echo $i
done
