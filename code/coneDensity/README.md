# cone density analysis
Clone the entire retinaTOMEAnalysis project. Issue the command `tbUseProject('retinaTOMEAnalysis')` to set up. 

These routines pick up where Rob Cooper's cone estimation process ends.

- processDensityMaps.m: Loads and filters the "merged" and "fovea" cone density maps that Rob's code produces. Saves some diagnostic images, and a `.mat` file with the cleaned density maps for each subject / eye / modality (merged and fovea).
