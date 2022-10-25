Overview
--------

The code in this folder conducts the empirical analysis and produces the results presented in ADD LINK TO PUBLISHED VERSION. The package consists of several scripts written in R and Python 3 that prepare datasets, apply statistical methods, and generate all figures and tables with results that are included in the manuscript and SI of the paper.

Data Availability
----------------------------

All data are publicly available at: ADD LINK TO ZENODO

Computational requirements
---------------------------

### Software Requirements

- Python 3.8.4
  - `numpy` 1.19.0
  - `pandas` 1.0.5
  - `scipy` 1.5.2
  - `statsmodels` 0.11.1
  - `sklearn` 0.23.1
  - `matplotlib` 3.2.2
  - `seaborn` 0.10.0
- R 3.6
  - `dplyr`
  - `lfe`
  - `MASS`
  - `glmnet`
  - `stargazer`
  - `boot`
  - `sandwich`

### Memory and Runtime Requirements

The code requires large memory because of the file of the dataset (approximately 100 GB).
Approximate time needed to reproduce the analyses on a HPC with 100 GB memory and 4 kernels: 72 hours.

Description of individual scripts
----------------------------

 - `01a_prepare_data_distlags.py`: This file prepares a dataset for the panel regression. It calculates the growth rate of synthetic cases (see paper for details) and the lags of all weather variables. The prepared data is written in csv format to the `data` subfolder.
 - `01b_prepare_data_distlags_firstdiff.py`: This file prepares a dataset for the panel regression. It calculates the growth rate of synthetic cases (see paper for details) and the lags of all weather variables. Furthermore, it takes first-differences of all variables. The prepared data is written in csv format to the `data` subfolder.
 - `01c_prepare_data_distlags_secdiff.py`: This file prepares a dataset for the panel regression. It calculates the growth rate of synthetic cases (see paper for details) and the lags of all weather variables. Furthermore, it takes second-differences of all variables. The prepared data is written in csv format to the `data` subfolder.
 - `01d_prepare_data_distlags_bins_firstdiff_L10.py`: This file prepares a dataset for the panel regression with bins for all weather variables. It calculates the growth rate of synthetic cases (see paper for details) and the lags of all weather variables. Furthermore, it takes first-differences of all variables. The prepared data is written in csv format to the `data` subfolder.
 - `01e_create_file_modelselection.py`: This file creates another file `modelselection_allvariables.csv` that contains the combinations of weather variables that are examined with the panel regression model.
 - `02a_smooth-models_distlags_L05.R`: This file estimates the panel regression model with 5 lags. The results of the estimation (coefficients and model statistics) are written in csv format to the `results` subfolder.
 - `02b_smooth-models_distlags_L10.R`: This file estimates the panel regression model with 10 lags. The results of the estimation (coefficients and model statistics) are written in csv format to the `results` subfolder.
  - `02c_smooth-models_distlags_L20.R`: This file estimates the panel regression model with 20 lags. The results of the estimation (coefficients and model statistics) are written in csv format to the `results` subfolder.
 - `02c_plot_figure_01.py`: This file produces Figure 1 in the paper using the results of the two previous scripts.
 - `03a_smooth-models_distlags_L10_modelselection.R`: This file reads in the file `modelselection_allvariables.csv` and estimates the panel regression model for all the different combinations of variables. The results of the estimation (coefficients and model statistics) are written in csv format to the `results` subfolder.
 - `03b_plot_variable_selection_figure_S12.py`: This file reads in the estimation results of the previous script and produces Figure 12 in the SI.
 - `04a_binned-models_distlags_L10.R`: This file estimates the panel regression model with 10 lags and bins of all weather variables. The results of the estimation (coefficients and model statistics) are written in csv format to the `results` subfolder.
 - `04b_plot_bins_figure_S14.py`: This file reads in the estimation results of the previous script and produces Figure 14 in the SI.
 - `05a_extract_dependent.py`: This file extracts the dependent variable (growth rate of synthetic cases) from the main dataset for subsequent analysis.
 - `05b_stationarity_figure_S11.py`: This file reads in the prepared data from the previous script and applies statistical tests of stationarity. It produces Figure 11 in the SI.
 - `05c_multicollinearity.py`: This file examines multicollinearity by calculating the VIF of different models with and without fixed effects and with and without first- or second-differencing. The results are stored in csv format in the `results` subfolder.
 - `05c_overfitting.py`: This file reads in estimation results from various models and produces Figure MISS in the SI that shows the model fit using the within-R2 and a metric from a 10-fold model cross-validation.
 - `05d_make_table_vif.py`: This file reads in estimation results and produces Table MISS in the SI that shows the results of the analysis of multicollinearity.
 - `06_plot_distlags_figure13.py`: This file reads in estimation results from the model with 20 lags and produces Figure 13 in the SI which shows the estimated coefficients for different time lags.
 - `07b_smooth-models_contemp_mobility.R`: This file estimates a panel regression model to examine the instantaneous effect of weather on mobility. The results of the estimation (coefficients and model statistics) are written in csv format to the `results` subfolder.
 - `07c_smooth-models_distlags_mortality.R`: This file estimates a panel regression model to examine the effect of weather on mortality. The results of the estimation (coefficients and model statistics) are written in csv format to the `results` subfolder.
 - `07d_smooth-models_distlags_mobility.R`: This file estimates a panel regression model to examine the effect of mobility on the disease growth. The results of the estimation (coefficients and model statistics) are written in csv format to the `results` subfolder.
 - `08_ftest.py`: This file conducts F tests based on the estimation results of different models (e.g. to test for the joint significance of the coefficients of the weather variables).

Instructions to Replicators
---------------------------

- Download the `panel_all.csv` file from the zenodo repository above into the folder `data_raw`.
- Run all scripts in the order indicated by the file names (i.e. `01`, `02`, `03`, ...). This can also be achieved with the Makefile in the repository (`make clean; make all`).
- Some of the scripts store processed data in the folder `data` or intermediate results in the folder `results`.
- Once all scripts have finished, all tables and figures can be found in the respective folders `tables` and `figures`.
