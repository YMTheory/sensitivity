#!/bin/bash

FILEDIR="/p/lustre2/lenardo1/sensitivity_output/Mar1_Xe137_Study_merged-v11_D024/"

for SCALE in 00.1x 00.3x 01.0x 03.0x 10.0x 100.0x
do

   for FILENAME in $(ls "$FILEDIR"*"$SCALE"*"h5")
   do 
      if [ "$SCALE" = "01.0x" ] 
      then
        CRITICAL_LAMBDA_FILE=$(ls "../Rn222Study/CriticalLambdaCurves/critical_lambda_merged-v11_D024_"*"$SCALE"*"txt")
      elif [ "$SCALE" = "00.3x" ]
      then
        CRITICAL_LAMBDA_FILE=$(ls "CriticalLambdaCurves/critical_lambda_xe137study_merged-v11_D024_00.1x.txt") 
      else
        CRITICAL_LAMBDA_FILE=$(ls "CriticalLambdaCurves/critical_lambda_xe137study_merged-v11_D024_"*"$SCALE"*"txt")
      fi
      echo $CRITICAL_LAMBDA_FILE
      python ../Reanalyze90PercentLimitWithNewCriticalLambda.py $FILENAME $CRITICAL_LAMBDA_FILE $FILEDIR
   done
        
echo $i
done
