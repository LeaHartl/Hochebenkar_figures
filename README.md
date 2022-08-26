# Hochebenkar 2022: figures

proc_data_stack.py is the main file and calls figures_stack.py and rasterPlotsStack.py

## proc_data_stack: 
- reads csv files of block positions and turns them into a dictionary of geodataframes for each line
- computes xy and xyz displacement and annual movement
- contains function to write point coordinates with velocities to shapefile
- calls various plotting functions in:

## rasterPlotsStack: 
contains code to plot 
- figures with subplots of hillshades, velocity vectors, vertical surf. elev. change
- rotation figure, strain rate figure
- overview figure

## figuresStack 
- flowline figures, single blocks, time series

## data folder 
- data/rasters: needs to exist for file paths in scripts to work. download stacked tifs to put in this subfolder at: https://doi.org/10.5281/zenodo.7010292
- shapefile of reference lines defining the block profiles.
- shapefile of central flow line
- shapefile (points) of block positions
- 2022_mittelProfile.csv: times series of profile mean velocities (m/a) up to 2021.
- mean_velocities_block_monitoring.txt: time series of mean DSM derived velocities at the location of the block profiles. 
For code related to the velocity vectors, see:
https://github.com/thomaszieher/HRG_reanalysis
