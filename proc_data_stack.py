#! /usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

import geopandas as gpd
from shapely.geometry import Point, LineString

# import files with code to make figures:
import figures_stack as fgs
import rasterPlotsStack as rp

# turn off setting with copy warning - CAREFUL!
pd.options.mode.chained_assignment = None  # default='warn'

# ---some helper functions----
# read csv files exported from martin's excel file.


def readFiles(fname, L):
    data = pd.read_csv(fname,  delimiter=';', header=1)
    data = data[['Date', 'Stone', 'North [m]', 'East [m]', 'Elevation [m]']]
    data.columns = ['Date', 'Stone', 'x', 'y', 'z']

    for c in ['x', 'y', 'z']:
        data[c] = data[c].astype(float)
    data['dxyz'] = np.nan
    data['dxyz_a'] = np.nan
    data['dxy'] = np.nan
    data['dz'] = np.nan
    data['Line'] = L
    data['Date'] = pd.to_datetime(data['Date'])

    for s in data.Stone.unique():
        d = data.loc[data['Stone'] == s]
        d['dz'] = d['z'].diff()
        data.loc[data['Stone'] == s, 'dz'] = d['dz'].values

        d['dxy'] = np.sqrt((d['x'].diff()**2 + d['y'].diff()**2))
        data.loc[data['Stone'] == s, 'dxy'] = d['dxy'].values

        d['dxyz'] = np.sqrt((d['x'].diff()**2 + d['y'].diff()**2 + d['z'].diff()**2))
        data.loc[data['Stone'] == s, 'dxyz'] = d['dxyz'].values

        # print(d['Date'].diff().astype(int))
        d['dxyz_a'] = (d['dxyz'] / d['Date'].diff().dt.days) * 365
        data.loc[data['Stone'] == s, 'dxyz_a'] = d['dxyz_a'].values

    geometry = [Point(xy) for xy in zip(data.y, data.x)]
    gdf = gpd.GeoDataFrame(data, crs="EPSG:31254", geometry=geometry)
    return(data, gdf)


# compute annual displacement
def annualavg(val, date1, date2):
    delta = datetime.strptime(date2, "%Y-%m-%d")-datetime.strptime(date1, "%Y-%m-%d")
    val2 = (val/abs(delta.days))*365
    return(val2)


# convert to decimal degrees
def convertDMS(DMS):
    deg = [i.split('°')[0] for i in DMS]
    minutes = [i.partition('°')[2].partition('′')[0] for i in DMS]
    seconds = [i.partition('′')[2].partition('′′')[0] for i in DMS]
    deg = np.array(deg, dtype=np.float32)
    minutes = np.array(minutes, dtype=np.float32)
    seconds = np.array(seconds, dtype=np.float32)

    dec = (deg.astype(float) + minutes.astype(float)/60 + seconds.astype(float)/(60*60))
    return(dec)


# make shapefile with all stone coordinates
def make_shp(Lines):
    #make one shapefile of stone coordinates
    L01 = Lines['L0']['gdf'].append(Lines['L1']['gdf'])
    L012 = L01.append(Lines['L2']['gdf'])
    L0123 = L012.append(Lines['L3']['gdf'])
    L01234 = L0123.append(Lines['L4']['gdf'])
    L01234.dropna(axis=0, subset=['Date'], inplace=True)
    L01234.Date = L01234.Date.astype(str)

    L01234.to_crs(epsg=32632, inplace=True)
    L01234.to_file('data/Steine_32632.shp')
    L01234.to_file('data/Steine_32632.geojson', driver='GeoJSON')

# ----------------------------


# load shapefile of RG outline
bg = gpd.read_file('data/HEK2018.shp')
# set crs to match stone coordinates
bg.to_crs(epsg=31254, inplace=True)

# initiatlize dictionary of line data
Lines = {
        'L0': {
            'fname': 'data/HEK_Bewegung_seit1996/Linie 0-Table 1.csv',
        },
        'L1': {
            'fname': 'data/HEK_Bewegung_seit1996/Linie 1-Table 1.csv',
                    },
        'L2': {
            'fname': 'data/HEK_Bewegung_seit1996/Linie 2-Table 1.csv',
                    },
        'L3': {
            'fname': 'data/HEK_Bewegung_seit1996/Linie 3-Table 1.csv',
        },
        'L4': {
            'fname': 'data/HEK_Bewegung_seit1996/Längs-Table 1.csv',
        }
        }


# initiatlize DF for mean velocities (only for cross profiles in the later years)
meanv = pd.DataFrame(index=np.arange(1997, 2022), columns=['L0', 'L1', 'L2', 'L3'])
maxv = pd.DataFrame(index=np.arange(1997, 2022), columns=['L0', 'L1', 'L2', 'L3'])


demyrs = [1953, 1971, 1977, 1990, 1997, 2006, 2009, 2010, 2011, 2017, 2018, 2019, 2020, 2021]
dispyrs = [1971, 1977, 1990, 1997, 2006, 2009, 2010, 2011, 2017, 2018, 2019, 2020, 2021]

# these are the raster files (DSM, displacement, elevation change). files are stacked tifs where
# each band/layer is a year, in chronological order.
demsHEK = 'data/raster/dsmstack.tif'
disp = 'data/raster/velstack.tif'
hs = 'data/raster/hsstack.tif'
dhHEK = 'data/raster/ddsmstack.tif'

# add gdf to Lines dictionary and populate DF of mean vel.
for L in Lines:
    dat, Lines[L]['gdf'] = readFiles(Lines[L]['fname'], L)
    dat['year'] = dat['Date'].dt.year
    Lines[L]['dat'] = dat

    for y in meanv.index:
        meanv.loc[y, L] = dat.loc[dat['year'] == y, 'dxyz_a'].mean()
        maxv.loc[y, L] = dat.loc[dat['year'] == y, 'dxyz_a'].max()

## uncomment look at values, print to csv if needed.
# print(meanv)
# print(maxv)
# maxv.to_csv('data/maxv.csv')


# deal with reference points: import coordinates from csv file, make a geodataframe
dat = pd.read_csv('data/fixpunkte.csv')
dat['decx'] = convertDMS(dat['startx'].astype(str).values)
dat['decy'] = convertDMS(dat['starty'].astype(str).values)

dat = dat[['p', 'decx', 'decy']]
gdf = gpd.GeoDataFrame(
    dat, geometry=gpd.points_from_xy(dat['decx'], dat['decy']))
gdf.set_crs(epsg=4326, inplace=True)
gdf.to_crs(epsg=31254, inplace=True)
gdf['y'] = gdf.geometry.y
gdf['x'] = gdf.geometry.x


## make another gdf and turn ref. points into lines, write to shapefile (adjust commenting out of lines as needed)
# LS0 = LineString([(gdf.loc[gdf['p'] == 'P0S', :].x, gdf.loc[gdf['p'] == 'P0S', :].y), (gdf.loc[gdf['p'] == 'P0E', :].x, gdf.loc[gdf['p'] == 'P0E', :].y)])
# LS1 = LineString([(gdf.loc[gdf['p'] == 'P1S', :].x, gdf.loc[gdf['p'] == 'P1S', :].y), (gdf.loc[gdf['p'] == 'P1E', :].x, gdf.loc[gdf['p'] == 'P1E', :].y)])
# LS2 = LineString([(gdf.loc[gdf['p'] == 'P2S', :].x, gdf.loc[gdf['p'] == 'P2S', :].y), (gdf.loc[gdf['p'] == 'P2E', :].x, gdf.loc[gdf['p'] == 'P2E', :].y)])
# LS3 = LineString([(gdf.loc[gdf['p'] == 'P3S', :].x, gdf.loc[gdf['p'] == 'P3S', :].y), (gdf.loc[gdf['p'] == 'P3E', :].x, gdf.loc[gdf['p'] == 'P3E', :].y)])

# d = {'col1': ['L0', 'L1', 'L2', 'L3'], 'geometry': [LS0, LS1, LS2, LS3]}
# df = gpd.GeoDataFrame(d, crs="EPSG:31254")
# df.to_file('data/ref_lines.shp')

## deal with reference lines, get points where they intersect the centerlines
# (import shapefile generated in the lines above).
ref_lines = gpd.read_file('data/ref_lines.shp')
ref_lines.to_crs(epsg=32632, inplace=True)
centerline = gpd.read_file('data/centerline_1.shp')

intersect = ref_lines.intersection(centerline.geometry.iloc[0])

d0 = centerline.geometry.iloc[0].project(intersect[0])
d1 = centerline.geometry.iloc[0].project(intersect[1])
d2 = centerline.geometry.iloc[0].project(intersect[2])
d3 = centerline.geometry.iloc[0].project(intersect[3])

# this contains the distance along the flowline for each profile
ds = [d0, d1, d2, d3]

# get the mean profile velocities for the entire block time series from csv file and format into nicer dataframe
meanV_all1 = pd.read_csv('data/2022_mittelProfile.csv', delimiter=';')
newix = np.arange(1953, 2022)
meanV_all = meanV_all1.copy()
meanV_all.set_index('year', inplace=True, drop=True)
meanV_all = meanV_all.reindex(newix, method='ffill')

# make percentage anomalies of velocity time series
allV = meanV_all.mean()

anomV = meanV_all.copy()
for c in anomV.columns:
    print(meanV_all[c].mean())
    anomV[c] = 100 * (meanV_all[c]) / meanV_all[c].mean()

#print(anomV.tail(30))

# ---------make shapefile of the block lines if needed------------
# make_shp(Lines)

# ---------plotting starts here-----------

# makes a few different (messy) timeseries plots of mean block profile velocities and BCF
# also returns data frame of Bdf per year along flowline, needed for time series plot
Bdf1, Bdf1_DSM = fgs.BCF_TS(meanV_all1, ds, demsHEK, demyrs)
## time series plot mean block profiles, dsm velocities at profile locations, BCF (fig 3)
fgs.timeseriesBoth4(meanV_all1, Bdf1, Bdf1_DSM, demyrs)
## main flowline plot with inset axis (fig 6)
fgs.FlowlineAll_inset(demsHEK, ds)
## Scarp plot (fig 8)
fgs.FlowlineScarp(demsHEK, ds)
## velocity along flowline from disp rasters, plot for supplement.
fgs.FlowlineVel(demsHEK, disp, ds, dispyrs)
# # subplots of single blocks (2016-2021) (fig 9)
fgs.SingleBlocks(Lines, gdf)
# makes flowline plot zoomed in on terminus for supplement.
fgs.FlowlineTerminus(demsHEK)

## raster plots:
## strain rates & map of blocks in P1 and P2 (fig 10)
# rp.figMap(Lines, bg, demsHEK)
## subplots of hillshades (fig 5)
# rp.plotHS(demsHEK, demyrs)
## subplots of velocity vectors (fig 4)
# rp.plotVectors2(demsHEK)
## subplots of elevation change (fig 7)
# rp.plotElevchange(demsHEK)
## overview plot (part of fig 1)
# rp.overview(Lines, demsHEK, ref_lines)
## rotation plot (pitch angles over hillshades and surf. elev. change, fig 11)
#rp.plotRotation(demsHEK, dhHEK)


plt.show()



