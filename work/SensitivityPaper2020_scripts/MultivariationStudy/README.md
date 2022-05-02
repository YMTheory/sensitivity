# Multi-parametric Sensitivity Study and Sensitivity-derived Requirements

These scripts expand the sensitivity studies for the 2021 Sensitivity paper to the DNN topological discriminator and account for varying multiple parameters simultaneously. 

A technical report with details and results is available on the project confluence. Some technical details for running the codes are provided here.

A copy of the results is stored on LLNL's HPSS at `/proj/nexo/Sensitivity2020/multivarstudy/`.

## Inputs

- Sensitivity 2020 Merged TTrees (on HPSS storage at `/proj/nexo/Sensitivity2020/merged-v11-mcid-labels/`)
- Materials Database D-024 v16

## Analysis steps

- First, use the `dnn_smoothing.ipynb` notebook to determine suitable smearing parameters for the DNN distribution, plot the resulting ROC curves, and compute the corresponding background misidentification rates and DNN cut values.
- Create histogram files (if not already available) for the necessary combinations of energy resolution and DNN response using `SubmitNewHistFromTTrees.py`. Do this using the config files for both standard and fine binning. 
- Created the Components Tables for all the histogram files from the previous step using `CreateComponentsTableFromHistograms.py`. This pulls from the database. Note histogram input and components output folders used here.
- Run the sensitivity calculations for the set of parameters of interest using `SubmitNew90CLSensitivityJobs_MultivarStudy_LLNL.py` to submit jobs running `Compute90PercentLimit_WilksApprox_MultivarStudy.py` to the cluster.
- Analysis of the results and plotting is performed in the `multivar analysis full.ipynb` notebook