# ETNP ODZ 50 year time series

"from https://zenodo.org/record/6519188#.YozZ0ZPMIdo" B.G. 2022-06-22 

This upload contains the datasets for Evans et al. (2022), "Natural variability and expansion of the nitrogen deficit within the Eastern Tropical North Pacific Oxygen Deficient Zone".
The primary data product here are 8 GLODAP-corrected cruises through the ETNP ODZ covering 50 years. This upload also includes the datasets and intermediate products for the integration
of fixed nitrogen loss as well as the CalCOFI data and Pescadero sediment core nitrogen isotope data. These additional data products are in subfolders for organization.

# ETNP ODZ time series metadata

The Excel sheet "ETNP ODZ file metadata.xlsx" contains information about the shorthand and Cruise IDs for each cruise. It also has the key for each data group for the integration calculations.
The second tab on this Excel sheet contains the GLODAP correction factors.

# Integration calculations folder

The Excel sheet "allCruises_updated_plus_PotDens_20220405.xlsx" contains the most updated data for this dataset. It identifies cruises based on "data group", which is less transparent, hence
why we include "All 8 Cruises_adjusted including KM1919 .csv" on the front of this upload. In addition to performing GLODAP corrections, this dataset set O2 = 0 when O2 < 10 μmol kg-1 and
NO2 > 0.1 μmol kg-1.

The file "cruiseSecQC_Integrations_final_intermediate_allCruises_updated20220405.mat" contains the data after GLODAP-correction.
The file "cruiseSecQC_Integrations_final_intermediate2_allCruises_updated20220405.mat" contains to the gridded data before integration.
The file "cruiseSecQC_Integrations_final_allCruises_updated20220405_Talias_integrations.mat" contains the integation results for the following integration subsets.
Method A: Density > 24.75 : 1000 meters
Method B: Density > 24.75 : Density < 27.2
Method C: Density > 24.75 : Density < 26.2
Method D: Density > 26.2 : Density < 26.8
Method E: Density > 26.8 : Density < 27.2
Method F: Density > 27.2 : 1000 meters

The file "cruiseSecQC_Integrations_final_allCruises_updated20220405.mat" contains the integration results for the older integration subsets.
Method A: Oxygen < 50 : 1000 meters
Method B: Oxygen < 50 : Density < 26.6
Method C: Density > 26.6 : 1000 meters
Method D: Density > 24.75 : 1000 meters
Method E: Density > 24.75 : Density < 26.6
Method F: Oxygen < 10 : 1000 meters
Method G: Oxygen < 10 : Density < 26.6
Method H: Oxygen < 50
Method I: Oxygen < 10

The "results_20220407.xlsx" file contains the integration results from both the older subsets and the Talia subsets as well as the GLODAP QC.

# eOMP folder

The eOMP folder contains the files that eOMP was used on and the outputs from this package. Each csv file was converted into a .mat file to read into the eOMP script. The file "All 8 Cruises_adjusted_omp.csv"
was converted to the "ETNP_cruises.mat" file, which was used for methods optimization. Each of the cruises was extracted as a csv file then converted to a mat file for eOMP calculation. 

The names in these columns were converted to text strings that the eOMP package could read. The key for these values is as follows:
oxy - oxygen/umol kg-1
ph - phosphate/umol kg-1
si - silicate/umol kg-1
ni - nitrate/umol kg-1
press = pressure/dbar
sal - absolute salinity/g kg-1
ptemp - conservative temperature/C
pdens - potential density/kg m-3
pvort - potential density/kg m-3

As for the outputs, variables in the output are saved in the same way. The matrix A contains the water mass distributions, in the following order: "13CW, NEPIW, AAIW, ESW, uPSUW".
The matrix "biogeo" contains the amount of remineralization calculated in each sample, and it's currently in units of carbon equivalents, as defined by 106*PO4 equivalents.

# CalCOFI and Pescadero data folder

Pescadero data, as downloaded from the SI of Tems et al. (2016), is saved in "Basin_data.mat". This file also contains the Santa Monica sediment core record from Tems et al. (2015), as shared by Will Berelson.
The file "calCOFI_6_var_mtrx_quarterly_dz05.nc" is a netcdf file containing the CalCOFI dataset, which Isaac Schroeder shared with me. In this dataset, data from 2020 and 2021 are not quality controlled and 
was not used in this analysis.

The following information describes the data structure of this netcdf. This is from Evans et al. (2020), focusing on the CCS, and this analysis only went to 2018, unlike this new dataset. 

% # -------------------------NetCDF Coordinates and Variables-----------------------
% # Dimensions:    (depth: 101, quarter: 4, station: 66, variable: 6, year: 70)
% # Coordinates:
% #   * station    (station) object '767_490' '767_510' ... '933_1100' '933_1200'
% #   * year       (year) int64 1949 1950 1951 1952 1953 ... 2015 2016 2017 2018
% #   * quarter    (quarter) int64 1 2 3 4
% #   * depth      (depth) int64 0 5 10 15 20 25 30 ... 470 475 480 485 490 495 500
% #   * variable   (variable) object 'R_POTEMP' 'R_SALINITY' ... 'R_SIO3' 'R_NO3'
% # Data variables:
% #     varQ_mtrx  (station, year, quarter, depth, variable) float64 ...
% #     monQ_mtrx  (station, year, quarter) float64 ...
% #     lon_vec    (station) float64 ...
% #     lat_vec    (station) float64 ...
% #     sttn       (station) float64 ...
% #     line       (station) float64 ...
% # ------------------------------------------------------------------------------

The file "CCS_data.mat" contains the netdcf decompressed, in a matrix named "varQ_mtrx", as well as the mean seasonal oxygen data, in the matrix "O2".

# Final notes

These datasets accompany the paper Evans et al. (2022), "Natural variability and expansion of the nitrogen deficit within the Eastern Tropical North Pacific Oxygen Deficient Zone".
There is a second upload on Zenodo containing the code for processing these datasets, and provides the scripts for converting between these data products. If you have questions, 
please reach out to Allan Devol or Natalya Evans.



