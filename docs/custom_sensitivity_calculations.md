# Running your own sensitivity calculations

Once you have [gone through the tutorial](https://github.com/nEXO-collaboration/sensitivity/blob/main/work/Sensitivity%20Code%20Tutorial.ipynb) 
and understand what the code is doing, the next step is to run your own large-scale calculation. 
There are a few key pieces to this puzzle:
- **The YAML config file**, which defines how the histograms are binned, what cuts are applied to the data when creating the PDFs, and how the components are grouped.
- **The ComponentsTable HDF5 file**, which contains the PDFs as well as a whole host of metadata that allows the code to appropriately scale the PDFs to create a background/signal model.
- **Calculation scripts**, which contain contain the actual code to run the calculations and store the output.

If you're running a brand new calculation, the steps are basically as follows:
1. Create your own YAML config file to set up the binning, cuts, and groupings that will be used to create the PDFs used in the fit
2. Create a file containing the histograms for each component using the [`BuildHistogramTableFromTrees.py`](https://github.com/nEXO-collaboration/sensitivity/blob/documentation_update/work/conversion_scripts/BuildHistogramTableFromTrees.py) script, then create a **ComponentsTable** using the [`CreateComponentsTableFromHistograms.py`](https://github.com/nEXO-collaboration/sensitivity/blob/documentation_update/work/conversion_scripts/CreateComponentsTableFromHistograms.py) script. 
3. Run the sensitivity calculation using a script like [`Compute90PercentLimit_WilksApprox_RnScaling_Example.py`](https://github.com/nEXO-collaboration/sensitivity/blob/documentation_update/work/example_scripts/Compute90PercentLimit_WilksApprox_RnScaling_Example.py)
4. Analyze the data to get distributions of the upper limits and compute the median sensitivity.

Details about each of these steps is provided below. 

<br/> 

## The YAML config file

The YAML file contains information that configures the fit. As an example, we can look at the configuration file used for the baseline sensitivity estimates in the 2020 publication: [```Sensitivity2020_Optimized_DNN_Standoff_Binning_version1.yaml```](https://github.com/nEXO-collaboration/sensitivity/blob/main/work/config/Sensitivity2020_Optimized_DNN_Standoff_Binning_version1.yaml)
Let's go section by section and break down what's in the file.

| **Section** | **Description** |
| ----------- | --------------- |
| `SummaryInputFile` and `RawPDFInputDir` | These sections are actually deprecated, and you can ignore them. I've retained them here for backwards compatibility with the 2017 analysis. (BL, July 2022) |
|  `FitAxes` | This section defines the "variable space" where we actually fit the data. In the 2020 sensitivity paper, this was a 3D space where "DNN", "Energy" and "Standoff" defined the three axes. Each axis needs to be explicitly defined here. You can see a few options here, that should be mostly self-explanatory. The two that are maybe not so obvious are: <ul><li> `Title`: this is where you, as a user, can choose the name of the variable as it will be referred to in the sensitivity code. </li><li> `VariableName`: this is the name of the variable *as it is defined in the ROOT data structure that is produced by the analysis framework*. This is required for building the histograms that then become the PDFs.</li></ul> <br/> The binning along each axis can either be `Linear` or `Custom`, the latter of which allows the user to manually set the edges of each bin. <br/><br/> **Note:** While the standard (as of July 2022) is to use the three variables DNN, Standoff, and Energy, there are some examples of fits with different numbers of variables in the repository; for instance, when we studied [how our sensitivity changes using only a 2D fit](https://github.com/nEXO-collaboration/sensitivity/blob/documentation_update/work/config/Sensitivity2020_Optimized_DNN_Standoff_Binning_version1_energy%2Bdnn.yaml). |
| `Cuts` | This section defines the cut that get applied to the data when constructing the PDFs. As of July 2022, there are three supported types of cuts: <ul><li> **Boolean**, which is just a yes/no (for instance, did the event pass some threshold) </li><li> **1D**, which is just what it sounds like (is variable `X` below or above a certain value for this event?) </li><li> **2DLinear**, which defines a linear cut along two different dimensions (we use this for e.g. the charge/light ratio cuts to remove skin events) </li></ul> |
| `CustomScalingFactors` and `CustomSpecificActivities` | Deprecated; ignore these. (BL, July 2022) |
| `Group Assignments` | This section assigns each of the **components** to a **group**, out of which we will build the PDFs for fitting. For instance, all the components that make up the Outer Cryostat are assigned to the "Far" group, while all the components that make up the <sup>238</sup>U background from the SiPM arrays are assigned to the "Internals_U238" group. <br><br> Components can be explicitly omitted from the analysis by assigning them to the "Off" group. For instance, in the 2020 sensitivity calculation we ignored any contributions from <sup>137</sup>Cs, so all of these components were "Off". <br/><br/> **Note:** *all* components in the background model need to be assigned to a group. If there is a component in the materials DB that does not get assigned a group in the YAML file, the code will crash when trying to create the ComponentsTable. |

<!---

#### `SummaryInputFile` and `RawPDFInputDir`:
These sections are actually deprecated, and you can ignore them. I've retained them here for backwards compatibility with the 2017 analysis. (BL, July 2022)

#### `FitAxes`:
This section defines the "variable space" where we actually fit the data. In the 2020 sensitivity paper, this was a 3D space where "DNN", "Energy" and "Standoff" defined the three axes. Each axis needs to be explicitly defined here. You can see a few options here, that should be mostly self-explanatory. The two that are maybe not so obvious are:
* `Title`: this is where you, as a user, can choose the name of the variable as it will be referred to in the sensitivity code. 
* `VariableName`: this is the name of the variable *as it is defined in the ROOT data structure that is produced by the analysis framework*. This is required for building the histograms that then become the PDFs.

The binning along each axis can either be `Linear` or `Custom`, the latter of which allows the user to manually set the edges of each bin.

> **Note:** While the standard (as of July 2022) is to use the three variables DNN, Standoff, and Energy, there are some examples of fits with different numbers of variables in the repository; for instance, when we studied [how our sensitivity changes using only a 2D fit](https://github.com/nEXO-collaboration/sensitivity/blob/documentation_update/work/config/Sensitivity2020_Optimized_DNN_Standoff_Binning_version1_energy%2Bdnn.yaml).


#### `Cuts`:
This section defines the cut that get applied to the data when constructing the PDFs. As of July 2022, there are three supported types of cuts:
* **Boolean**, which is just a yes/no (for instance, did the event pass some threshold)
* **1D**, which is just what it sounds like (is variable `X` below or above a certain value for this event?)
* **2DLinear**, which defines a linear cut along two different dimensions (we use this for e.g. the charge/light ratio cuts to remove skin events)



#### `CustomScalingFactors` and `CustomSpecificActivities`:
Deprecated; ignore these. (BL, July 2022)


#### `Group Assignments`:
This section assigns each of the **components** to a **group**, out of which we will build the PDFs for fitting. For instance, all the components that make up the Outer Cryostat are assigned to the "Far" group, while all the components that make up the <sup>238</sup>U background from the SiPM arrays are assigned to the "Internals_U238" group. 

Components can be explicitly omitted from the analysis by assigning them to the "Off" group. For instance, in the 2020 sensitivity calculation we ignored any contributions from <sup>137</sup>Cs, so all of these components were "Off".

> **Note:** *all* components in the background model need to be assigned to a group. If there is a component in the materials DB that does not get assigned a group in the YAML file, the code will crash when trying to create the ComponentsTable.

--->

<br/>

## The **ComponentsTable** HDF5 file

The **ComponentsTable** is a `pandas.DataFrame` object stored in an HDF5 file format. It can be opened using the `pandas.from_hdf()` function. As explained in the tutorial, the **ComponentsTable** combines three pieces of information:
1. Data from the MaterialDB, such as the mass, material, and specific activity of each component
2. The event distribution, stored as a binned **histogram**
3. The **group**, as defined in the YAML config file

> As an example, the ComponentsTable used in the 2020 baseline sensitivity estimate is stored on the LLNL machines here: ```/p/vast1/nexo/sensitivity2020/pdfs/component_tables/ComponentsTable_D-024_merged-v11_Optimized_DNN_Standoff_Binning_version1.h5 ```


Before creating a new **ComponentsTable**, you need to generate a histogram file. There are two scripts, located in the repository at `sensitivity/work/conversion_scripts`:
-  [**`BuildHistogramTableFromTrees.py`**](https://github.com/nEXO-collaboration/sensitivity/blob/documentation_update/work/conversion_scripts/BuildHistogramTableFromTrees.py), which converts the merged ROOT trees from the data processing/reconstruction into python-readable histograms from which we can build the PDFs. *This is the step where the code reads the binning and axis definitions from the YAML file.* You need to generate histograms *before* generating a **ComponentsTable**.
- [**`CreateComponentsTableFromHistograms.py`**](https://github.com/nEXO-collaboration/sensitivity/blob/documentation_update/work/conversion_scripts/CreateComponentsTableFromHistograms.py), which combines the **histograms** from the file with information from the MaterialDB to create a **ComponentsTable**. 

For the 2020 sensitivity estimate, the processed/reconstructed ROOT data can be found at:
```/p/vast1/nexo/merged-v11-mcid-labels/```

> **Note:** Generating the histograms can take up to 45 minutes, especially when the binning is chosen pretty finely. It will take the most time on PDFs like the 0nu and 2nu, for which we have lots of Monte Carlo events that pass our cuts. Creating the **ComponentsTable**, on the other hand, should only take a minute or so.

<br/> 

## Run scripts and output data structure

We've included some example scripts in the directory [`sensitivity/work/example_scripts`](https://github.com/nEXO-collaboration/sensitivity/tree/documentation_update/work/example_scripts) to get you started. 

For now, we'll focus on [`Compute90PercentLimit_WilksApprox_RnScaling_Example.py`](https://github.com/nEXO-collaboration/sensitivity/blob/documentation_update/work/example_scripts/Compute90PercentLimit_WilksApprox_RnScaling_Example.py).

> **Note:** We've also included an example script for submitting large batches of jobs to the LLNL cluster, [`SubmitNew90CLSensitivityJobs_LLNL_RnScaling_Example.py`](https://github.com/nEXO-collaboration/sensitivity/blob/documentation_update/work/example_scripts/SubmitNew90CLSensitivityJobs_LLNL_RnScaling_Example.py)


In this example, the run script is intended to compute nEXO's sensitivity with the option to scale up/down the contribution of <sup>222</sup>Rn-induced backgrounds. The anatomy of this script is as follows:
- **Lines 5-22:** Define the arguments:
    - `job_id_num` is just an identifier that gets added to the output filename; useful when running lots of jobs in parallel
    - `num_datasets` is the number of toy datasets that the code will run
    - `input_table` is the **ComponentsTable** that contains the PDFs and associated data
    - `output_dir` is where the output files will go
    - `rn222_scale_factor` is a user-defined number that will scale the <sup>222</sup>Rn background PDF up or down
- **Lines 24-65:** Define a function to grab the 90% CL upper limit using a quadratic interpolation between diffent hypotheses; this is the same method used in the [tutorial notebook](https://github.com/nEXO-collaboration/sensitivity/blob/documentation_update/work/Sensitivity%20Code%20Tutorial.ipynb)
- **Lines 69-172:** Set up the workspace and input settings; print out background breakdown, etc.
- **Lines 173-190:** Set up calculation loop
- **Line 191:** Start loop over toy datasets
- **Lines 193-246:** Rebuild model for each toy dataset (required for fluctuating the radioassay numbers), and apply constraints
- **Lines 249-280:** Loop over different signal hypotheses, perform fits, and calculate the profile likelihood ratios
- **Lines 283-333:** Store output from the fitter in a data structure that can be used to perform analysis later
- **Lines 337-347:** Write data to output file.

The data that gets saved to the output file can be modified by the user, but below we list the variables saved in the 2020 baseline calculation of our 90% CL sensitivity. **Note: One row of data is saved *per toy dataset* **


#### Output data for 2020 90% CL sensitivity analysis
| **Variable** | **Description** |
|:------------| :---------------|
| `num_signal` | array containing the values of the signal hypothesis for various hypothesis tests |
| `lambda` | array containing the value of the log profile likelihood ratio for each hypothesis test |
| `fixed_fit_converged` | array of booleans indicating whether or not the fit (with fixed signal) converged for each hypothesis test |
| `fixed_fit_acc_covar` | array of booleans indicating whether or not the fit (with fixed signal) produced a good covariance matrix for each hypothesis test |
| `fixed_fit_parameters` | values of free parameters for the hypothesis test fit (with fixed signal) |
| `fixed_fit_errors` | uncertainties on free parameters for the hypothesis test fit (with fixed signal) |
| `num_iterations` | array of number of tries it took to get a fit (with fixed signal) which converged properly... maximum is 10, after which it stops trying |
| `best_fit_converged` | boolean indicating whether or not the fit (with floating signal) converged |
| `best_fit_covar` | boolean indicating whether or not the fit (with floating signal) produced a good covariance matrix |
| `best_fit_iterations` | number of tries it took to get a fit (with floating signal) which converged properly... maximum is 10, after which it stops trying |
| `best_fit_parameters` | values of free parameters after the fit (with floating signal) |
| `best_fit_errors` | uncertainties on free parameters after the fit (with floating signal) |
| `best_fit_nll` | value of negative log likelihood for the fit (with floating signal) |
| `input_parameters` | initial values of the free parameters for the fit (with floating signal) |
| `90CL_crossing` | the signal hypothesis where the profile-likelihood-ratio test statistic crosses the 90% CL threshold, assuming Wilks' theorem |


<br/>

## Analyzing the output to compute actual sensitivity and discovery potentials

As an example of how to analyze the output data, please refer to some example notebooks from the 2020 sensitivity evaluation, for instance our analysis of how our sensitivity scales with Rn222 background:
* [work/SensitivityPaper2020_scripts/Rn222Study/Rn222 Analysis.ipynb](https://github.com/nEXO-collaboration/sensitivity/blob/documentation_update/work/SensitivityPaper2020_scripts/Rn222Study/Rn222%20Analysis.ipynb)




