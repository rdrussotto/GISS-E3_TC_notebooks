# GISS-E3_TC_notebooks

Jupyter notebooks, Python scripts, and data files used to make figures analyzing tropical cyclones (TCs) in the NASA GISS-E3 global climate model.
Much of this code was adapted from Matlab code written by Jeffrey Strong.

This code is for a paper to be submitted to Journal of Advances in Modeling Earth Systems in 2021 by 
Rick Russotto, Jeffrey Strong, Suzana Camargo, Adam Sobel, Gregory Elsaesser, Maxwell Kelley, Anthony Del Genio, Yumin Moon, and Daehyun Kim. 
The code and data will be uploaded to Zenodo upon acceptance of the paper. 

Vesions of libraries used for paper: Python 3.6.12, NumPy 1.19.2, Matplotlib 3.2.0 (including Basemap), Pandas 1.1.5, XArray 0.16.2, SciPy 1.5.2

### Notebooks used to make figures:

Figure 1: Plot_Tracks.ipynb

Figure 2: Plot_Density.ipynb

Figures 3-7: Plot_Statistics.ipynb

Figures 8-9: Plot_Storm_Tangential.ipynb

Figure 10: Plot_Radial_Profiles.ipynb

Figures 11-12: Sensitivity_Test_Plots.ipynb

Figure 13: Intermediate_Error_Plots_v2.ipynb

Figures 14, S2: Intermediate_Error_Plots_v3.ipynb

### Additional dependencies:

#### TC track preprocessing scripts: 
Read_Zhao_TCs.ipynb (for modeled TCs)

Preprocess_IBTrACS.ipynb (for observed TCs)

#### Storm-centered diagnostics:
Storm_Centered_Plots.py

Storm_Centered_V1_Prec.py

Retrieve_Utan_Vtan.py

#### Radial-height profiles:
Figures 10 and S1 depend on Matlab scripts (not included here) developed by Yumin Moon for paper in Journal of Climate, Moon et al. (2020), https://doi.org/10.1175/JCLI-D-19-0172.1, and adapted by Jeffrey Strong for this study.

#### Sensitivity test analysis:
Figures 11-14 and S2 depend on Matlab scripts developed by Jeffrey Strong (not included here) for processing climate variables in the model results and comparing to observations.

### Data files:

#### Radial profiles:
YuminEtAl_Figs.mat

#### Sensitivity test line plots (Figures 11-12): 
JS_C180_ParamSensTest.mat

JS_C180v2_ParamSensTest.mat

#### Values plotted in color matrix plots for Figures 13, 14, S2:
(Post NetCDF files here after saving)

### Data files not included due to file size (but will be in final repository):

#### Modeled TC tracks:
zhao_tracks_v1.nc

zhao_tracks_v2.nc 

(Observed TC data from IBTrACS project available at https://www.ncdc.noaa.gov/ibtracs/index.php?name=ib-v4-access)

#### Storm-centered TC quantities: 
storm_centered_peak_intensity_v1.nc

storm_centered_peak_intensity_v2.nc

storm_centered_reference_winds_v1.nc

storm_centered_reference_winds_v2.nc

storm_centered_peak_intensity_prec_v1.nc
