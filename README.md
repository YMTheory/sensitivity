# sensitivity

This repository contains the code used to calculate the sensitivity of nEXO. It performs the following functions:
1. Convert the ROOT files from the G4/reconstruction/merging step into python-readable histograms.
2. Download radioassay data from the MaterialDB, and associate it with the correct histograms.
3. Combine the simulated backgrounds from each component into grouped PDFs.
4. Generate toy datasets and perform a profile likelihood ratio analysis to estimate our sensitivity and/or discovery potential.

# Dependencies
The code requires the following external libraries
* [uproot](https://github.com/scikit-hep/uproot5) -- Main branch is currently compatible with uproot4.
* [pandas](https://pandas.pydata.org/) - Tested up to v0.24
* [xlrd](https://xlrd.readthedocs.io/en/latest/) - Tested up to v1.2
* [iminuit](https://pypi.org/project/iminuit/) - Tested up to v2.8.4
* [histlite](https://histlite.readthedocs.io/en/latest/)
* [cloudant](https://github.com/cloudant/python-cloudant)
* [pyyaml](https://pyyaml.org/)
* tables (if you are on a recent mac, see https://github.com/PyTables/PyTables/issues/828)  

All of the above can be installed using `pip install`.

If you are running on LLNL, I recommend creating a virtual environment by following the instructions under "Installing Your Own Site Packages" here: [https://hpc.llnl.gov/software/development-environment-software/python](https://hpc.llnl.gov/software/development-environment-software/python)

# Getting started

### Low-level tutorial

If you are new to the sensitivity calculation, we recommend you start with the tutorial. 
This can be found in [https://github.com/nEXO-collaboration/sensitivity/blob/main/work/Sensitivity%20Code%20Tutorial.ipynb](`sensitivity/work/Sensitivity Code Tutorial.ipynb`).
The tutorial is designed to walk through how the sensitivity calculation works, step-by-step. 

### Running your own sensitivity calculations

If you're familiar with how the calculation works, you'll want to run large ensembles of fits to toy datasets. There is some documentation to get you started here:

- [**Running your own sensitivity calculations**](https://github.com/nEXO-collaboration/sensitivity/blob/documentation_update/docs/custom_sensitivity_calculations.md)


# About this repository
Mostly maintained by Brian (blenardo@stanford.edu) and Samuele (sangiorgio1@llnl.gov), as of Dec 30, 2021.
