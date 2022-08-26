# Hochebenkar 2022: figures

proc_data_stack.py is the main file and calls figures_stack.py and rasterPlotsStack.py

## proc_data_stack: 
- reads csv files generated from Martin's excel and turns them into a dictionary of geodataframes for each line
- computes xy and xyz displacement and annual movement
- contains function to write point coordinates with velocities to shapefile
- contains funtion to convert degrees/min/sec to decimal degrees
- sets up dictionary with file names of DEMs (key = year, fname = filename)

## rasterPlotsStack: 
contains code to plot 
- figures with subplots of hillshades, velocity vectors, vertical surf. elev. change
- rotation figure, strain rate figure
- overview figure

## figuresStack 
- flowline figures, single blocks, time series

## data folder (needs update)
- subfolder DEMs: aligned, sometimes resampled versions of the various DEMs for use in plotting and other analysis. (subtracted GEOID from some - I think the CRS all match but would be good to double check if anything is shifted)
- Shapeflie of rock glacier outline, roughly mapped by me based on 2017/18 BEV DEM
- Shapefile of central flowline, generated with CenterLine python package and adjusted a little in QGIS
- fixpunkte.csv: coordinates (in DMS format) of fix points for the 4 cross profiles
- daten_christoph.csv: mean & max horizontal displacement values, copied into file from Klug et al., 2012.
- 2022_mittelProfile.csv: times series of profile mean velocities (m/a) up to 2021.
