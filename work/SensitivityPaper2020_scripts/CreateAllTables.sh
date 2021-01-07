#!/bin/bash

version=merged-v10b
LOCATION=/p/vast1/nexo/data/$version-mcid-labels/
COMMON=/p/vast1/nexo/sensitivity2020/pdfs/

Tags=(
  D-023
  D-024
)

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

for Config in "${ConfigFiles[@]}"; do
  echo $Config
  python BuildHistogramTableFromTrees.py ../config/Sensitivity2020_${Config}.yaml ${version}_${Config} $LOCATION ./
  for tag in "${Tags[@]}"; do
    python CreateComponentsTableFromHistograms.py ../config/Sensitivity2020_${Config}.yaml ${version}_${Config} Baseline2019_Histograms_${version}_${Config}.h5 ./ ${tag}
    mv ComponentsTable_${tag}_${version}_${Config}.h5 ${COMMON}/component_tables/
    chmod ug+rw ${COMMON}/component_tables/ComponentsTable_${tag}_${version}_${Config}.h5
    chgrp nexo ${COMMON}/component_tables/ComponentsTable_${tag}_${version}_${Config}.h5
  done
  mv Baseline2019_Histograms_${version}_${Config}.h5 ${COMMON}/histogram_files/
  chmod ug+rw ${COMMON}/histogram_files/Baseline2019_Histograms_${version}_${Config}.h5
  chgrp nexo ${COMMON}/histogram_files/Baseline2019_Histograms_${version}_${Config}.h5
done

