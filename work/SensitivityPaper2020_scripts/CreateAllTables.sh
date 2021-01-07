#!/bin/bash

# Author: Ako Jamil (ako.jamil@yale.edu)
# This script is intended to generate all needed sets of PDFs 
# for the sensitivity and discovery potential calculations given 
# any changes in the reconstruction output files 
# or binning or grouping of the PDFs defined in the Yaml config files  

# Version of the merged trees
version=merged-v10b

# Location of the ROOT files 
LOCATION=/p/vast1/nexo/data/$version-mcid-labels/

# Location where the PDFs will be saved
COMMON=/p/vast1/nexo/sensitivity2020/pdfs/

# Geometry tags corresponding to a model with Aurubis Copper (D-023)
# or electro-formed Copper (D-024) for all copper in the detector
Tags=(
  D-023
  D-024
)

# List of different Yaml config files to use for generating different sets of PDFs
ConfigFiles=(
  Binary_DNN_085
  Binary_DNN_0915
  BkgIndex_Binning
  Optimized_DNN_Standoff_Binning_version1_energy+dnn
  Optimized_DNN_Standoff_Binning_version1_energy+recon_standoff
  Optimized_DNN_Standoff_Binning_version1_energy_only
  Optimized_DNN_Standoff_Binning_version1_energy_only_with_binary_dnn_085
  Optimized_DNN_Standoff_Binning_version1
  Optimized_DNN_Standoff_Binning_version2
  Optimized_DNN_Standoff_Binning_version3
  Optimized_DNN_Standoff_Binning_version4
)

# Running the code
for Config in "${ConfigFiles[@]}"; do
  echo $Config

  # First create the corresponding histogram files from the merged trees, this usually takes a long time
  python BuildHistogramTableFromTrees.py ../config/Sensitivity2020_${Config}.yaml ${version}_${Config} $LOCATION ./

  # For each geometry tag create a corresponding components table, move it to the common location and change group access
  for tag in "${Tags[@]}"; do
    python CreateComponentsTableFromHistograms.py ../config/Sensitivity2020_${Config}.yaml ${version}_${Config} Baseline2019_Histograms_${version}_${Config}.h5 ./ ${tag}
    mv ComponentsTable_${tag}_${version}_${Config}.h5 ${COMMON}/component_tables/
    chmod ug+rw ${COMMON}/component_tables/ComponentsTable_${tag}_${version}_${Config}.h5
    chgrp nexo ${COMMON}/component_tables/ComponentsTable_${tag}_${version}_${Config}.h5
  done

  # Finally, move the histogram files to the common location and change group access
  mv Baseline2019_Histograms_${version}_${Config}.h5 ${COMMON}/histogram_files/
  chmod ug+rw ${COMMON}/histogram_files/Baseline2019_Histograms_${version}_${Config}.h5
  chgrp nexo ${COMMON}/histogram_files/Baseline2019_Histograms_${version}_${Config}.h5
done

