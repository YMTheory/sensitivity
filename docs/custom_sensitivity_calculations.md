# Running your own sensitivity calculations

Once you have [gone through the tutorial](https://github.com/nEXO-collaboration/sensitivity/blob/main/work/Sensitivity%20Code%20Tutorial.ipynb) 
and understand what the code is doing, the next step is to run your own large-scale calculation. 
There are a few key pieces to this puzzle:
1. **The YAML config file**, which defines how the histograms are binned, what cuts are applied to the data when creating the PDFs, and how the components are grouped.
2. **The ComponentsTable HDF5 file**, which contains the PDFs as well as a whole host of metadata that allows the code to appropriately scale the PDFs to create a background/signal model.
3. **Calculation scripts**, which contain contain the actual code to run the calculations and store the output.

We'll discuss each of these in some detail below. 


## The YAML config file

The YAML file contains information that configures the fit. As an example, we can look at the configuration file used for the baseline sensitivity estimates in the 2020 publication: [```Sensitivity2020_Optimized_DNN_Standoff_Binning_version1.yaml```](https://github.com/nEXO-collaboration/sensitivity/blob/main/work/config/Sensitivity2020_Optimized_DNN_Standoff_Binning_version1.yaml)
Let's go section by section and break down what's in the file.

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




