#! /usr/bin/env python3
import numpy as np
import pandas as pd
# import os
import matplotlib.pyplot as plt
# import fiona
# import matplotlib.colors as colors
# from matplotlib import cm
from matplotlib_scalebar.scalebar import ScaleBar
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as colors
import matplotlib.patheffects as path_effects
# import matplotlib.cbook as cbook
from matplotlib import cm
import geopandas as gpd
import rasterio
import earthpy.spatial as es
import richdem as rd
import xarray as xr


def getdata(f, i):
    with rasterio.open(f) as src:
        el = src.read(i+1)
        print(i+1)
        BBox = (src.bounds[0],  src.bounds[2], src.bounds[1],  src.bounds[3])
        # Set masked values to np.nan
        el[el < 0] = np.nan
        hs = es.hillshade(el, altitude=10, azimuth=30)

        gt = src.transform
        #PIXELSIZE
        piX = gt[0]
        piY = gt[4]
        #print(piX, piY)

    return(el, hs, BBox)


def getdata2(f, i):
    with rasterio.open(f) as src:
        el = src.read(i+1)
        BBox = (src.bounds[0],  src.bounds[2], src.bounds[1],  src.bounds[3])
        # Set masked values to np.nan
        # el[el < 0] = np.nan
        # hs = es.hillshade(el, altitude=10)

        gt = src.transform
        #PIXELSIZE
        piX = gt[0]
        piY = gt[4]
        #print(piX, piY)

    return(el, BBox)


def plotHS(demsHEK, demyrs):

    al = 0.9
    # fig, ax = plt.subplots(3, 5, figsize=(12, 6), sharex=True, sharey=True)
    fig, ax = plt.subplots(5, 3, figsize=(7, 10), sharex=True, sharey=True)
    ax = ax.flatten()
    # Open the DEM with Rasterio, make and plot hillshades with earthpy
    for i, d in enumerate(demyrs):
        #print(d)
        el, hs, BBox = getdata(demsHEK, i)
        ax[i].imshow(hs, cmap="Greys", alpha=al, extent=BBox)
        ax[i].grid('both')

    ax[0].set_xlim([51000, 51800])
    ax[0].set_ylim([188600, 189200])
    ax[10].set_xticks([51200, 51400, 51600])
    ax[10].set_yticks([188700, 188900, 189100])
    ax[-1].axis('off')

    for i, a in enumerate(ax[:-1]):
        lbl1 = ax[i].annotate('A', xy=(51200, 188770),  xycoords='data',
                       fontsize=15, horizontalalignment='left', verticalalignment='top',)
        lbl2 = ax[i].annotate('B', xy=(51160, 188900),  xycoords='data',
                       fontsize=15, horizontalalignment='left', verticalalignment='top',)
        lbl3 = ax[i].annotate(str(demyrs[i]), xy=(51750, 189000),  xycoords='data',
                       fontsize=15, horizontalalignment='right', verticalalignment='bottom',)
        lbl1.set_path_effects([path_effects.Stroke(linewidth=1.5, foreground='white'), path_effects.Normal()])
        lbl2.set_path_effects([path_effects.Stroke(linewidth=1.5, foreground='white'), path_effects.Normal()])
        lbl3.set_path_effects([path_effects.Stroke(linewidth=1.5, foreground='white'), path_effects.Normal()])

        scalebar = ScaleBar(1)# 1 pixel = 1 meter
        ax[i].add_artist(scalebar)

    # fig.subplots_adjust(hspace = 0, wspace = 0, bottom = 0.04, top=0.99, right =0.99 , left = 0.06)
    fig.subplots_adjust(hspace=0, wspace=0, bottom=0.04, top=0.99, right=0.99, left=0.10)
    fig.savefig('figs/hillshades.png', dpi=300)


def plotVectors2(demsHEK):
    al = 0.9
    # fig, ax = plt.subplots(3, 5, figsize=(12, 8), sharex=True, sharey=True)
    fig, ax = plt.subplots(4, 4, figsize=(8, 9), sharex=True, sharey=True)

    ax = ax.flatten()
    divider = make_axes_locatable(ax[-2])
    # set DEMs to be plotted
    demyrs = [1971, 1977, 1990, 1997, 2006, 2009, 2010, 2011, 2017, 2018, 2019, 2020, 2021]
    ttlyrs = ['1953-71', '1971-77', '1977-90', '1990-97', '1997-06', '2006-09', '2009-10', '2010-11',
              '2011-17', '2017-18', '2018-19', '2019-20', '2020-21']

    # Open the DEM with Rasterio, make and plot hillshades with earthpy
    for i, d in enumerate(demyrs):
        print(d)
        el, hs, BBox = getdata(demsHEK, i+1)
        im = ax[i].imshow(hs, cmap="Greys", alpha=al, extent=BBox)
        ax[i].grid('both')

    # add time step as annotation
    for i, t in enumerate(ttlyrs):
        tt = ax[i].annotate(t, xy=(51600, 188150),  xycoords='data',
                            fontsize=15, horizontalalignment='right', verticalalignment='top',)
        tt.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='black'))


    # set axis limits, ticks, etc
    ax[0].set_xlim([51000, 52000])
    ax[0].set_ylim([188000, 189150])
    ax[-1].axis('off')
    ax[-2].axis('off')
    ax[-3].axis('off')

    ax[9].set_xticks([51500])

    scalebar = ScaleBar(1) # 1 pixel = 1 meter
    ax[9].add_artist(scalebar)

    # get shapefiles with velocity vectors
    velvec = {1971: 'vel_vect_53-71.shp',
              1977: 'vel_vect_71-77.shp',
              1990: 'vel_vect_77-90.shp',
              1997: 'vel_vect_90-97.shp',
              2006: 'vel_vect_97-06.shp',
              2009: 'vel_vect_06-09.shp',
              2010: 'vel_vect_09-10.shp',
              2011: 'vel_vect_10-11.shp',
              2017: 'vel_vect_11-17.shp',
              2018: 'vel_vect_17-18.shp',
              2019: 'vel_vect_18-19.shp',
              2020: 'vel_vect_19-20.shp',
              2021: 'vel_vect_20-21.shp'}
    velfolder = '/Users/leahartl/Desktop/HEK/FINAL_RESULTS_GKWest/velocity_vectors/'

    # set color map for velocities
    cmap = cm.get_cmap('viridis_r', 15)

    norm = colors.Normalize(vmin=0, vmax=30)
    cbar_ax = fig.add_axes([0.5, 0.15, 0.4, 0.02])

    # get profile reference lines
    reflines = gpd.read_file('data/ref_lines.shp')
    # print(reflines.crs)

    # plot vectors 
    fts = 10
    for j, v in enumerate(velvec):
        vv = gpd.read_file(velfolder+velvec[v])
        v = vv.plot(ax=ax[j], alpha=0.9, column='velocity', vmin=0, vmax=30, cmap=cmap)
    # plot ref lines, annotate ref lines
    for j, d in enumerate(demyrs): 
        if d != 2018:   
            rf = reflines.plot(ax=ax[j], colors='k', linewidth=0.5)
            ax[j].annotate('P0', xy=(51500, 189110),  xycoords='data',
                       fontsize=fts,
                       horizontalalignment='right', verticalalignment='top',)
            ax[j].annotate('P1', xy=(51600, 188980),  xycoords='data',
                       fontsize=fts,
                       horizontalalignment='right', verticalalignment='top',)
            ax[j].annotate('P2', xy=(51720, 188760),  xycoords='data',
                       fontsize=fts,
                       horizontalalignment='right', verticalalignment='top',)
            ax[j].annotate('P3', xy=(51910, 188640),  xycoords='data',
                       fontsize=fts,
                       horizontalalignment='right', verticalalignment='top',)
        if d == 2018:
            rf01 = reflines.loc[(reflines.col1 == 'L0') | (reflines.col1 == 'L1')]
            rf01.plot(ax=ax[j], colors='k', linewidth=0.5)
            ax[j].annotate('P0', xy=(51500, 189110),  xycoords='data',
                       fontsize=fts,
                       horizontalalignment='right', verticalalignment='top',)
            ax[j].annotate('P1', xy=(51600, 188980),  xycoords='data',
                       fontsize=fts,
                       horizontalalignment='right', verticalalignment='top',)


    # add colorbar
    fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), cax=cbar_ax, orientation='horizontal')

    # adjust whitespace
    fig.subplots_adjust(hspace=0, wspace=0, bottom=0.04, top=0.99, right=0.99, left=0.1)
    # fig.subplots_adjust(hspace = 0, wspace = 0, bottom = 0.04, top=0.99, right =0.99 , left = 0.06)

    cbar_ax.set_xlabel('velocity (m/a)')

    fig.savefig('figs/velvectors2.png', dpi=300)


def plotElevchange(demsHEK):
    #fig, ax = plt.subplots(3, 5, figsize=(12, 8), sharex=True, sharey=True)
    fig, ax = plt.subplots(4, 4, figsize=(8, 9), sharex=True, sharey=True)

    ax = ax.flatten()
    divider = make_axes_locatable(ax[-2])

    # set DEMs to be plotted
    demyrs = [1971, 1977, 1990, 1997, 2006, 2009, 2010, 2011, 2017, 2018, 2019, 2020, 2021]
    ttlyrs = ['1953-71', '1971-77', '1977-90', '1990-97', '1997-06', '2006-09', '2009-10', '2010-11',
              '2011-17', '2017-18', '2018-19', '2019-20', '2020-21']

    # Open the DEM with Rasterio, make and plot hillshades with earthpy
    for i, d in enumerate(demyrs):
        print(d)
        el, hs, BBox = getdata(demsHEK, i+1)
        im = ax[i].imshow(hs, cmap="Greys", alpha=1, extent=BBox)
        ax[i].grid('both')

    for i, t in enumerate(ttlyrs):
        tt = ax[i].annotate(t, xy=(51600, 188150),  xycoords='data',
                            fontsize=15, horizontalalignment='right', verticalalignment='top',)
        tt.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='black'))

    ax[0].set_xlim([51000, 52000])
    ax[0].set_ylim([188000, 189150])
    ax[-1].axis('off')
    ax[-2].axis('off')
    ax[-3].axis('off')
    ax[11].set_xticks([51500])

    scalebar = ScaleBar(1)# 1 pixel = 1 meter
    ax[9].add_artist(scalebar)

    cbar_ax = fig.add_axes([0.5, 0.15, 0.4, 0.02])

    clrs = cm.get_cmap('RdBu', 12)

    norm = colors.Normalize(vmin=-2, vmax=2)
    alphas = [1, 1, 1, 1, 1, 0.2, 1, 1, 1, 1, 1, 1]
    dhHEK = 'data/raster/ddsmstack.tif'
    for i, d in enumerate(demyrs):
        el, BBox = getdata2(dhHEK, i)

        # this avoids plotting no data values.
        el[el < -20] = np.nan
        el[el > 20] = np.nan

        # adjust transparency for low change values so that hillshades are visible
        alphas = np.ones(el.shape)
        alphas[(el < 0.2) & (el > -0.2)] = 0.2
        
        im2 = ax[i].imshow(el, cmap=clrs, alpha=alphas, extent=BBox, vmin=-2, vmax=2)
        # ax[i].set_title(ttlyrs[i])
        ax[i].grid('both')
    cb = fig.colorbar(im2, cax=cbar_ax, orientation='horizontal', label='elevation change (m/a)')
    cbar_ax.set_xticks([-2, -1, 0, 1, 2])

    cbar_ax.set_xticklabels(['< -2', '-1', '0', '1', '> 2'])

    # get profile reference lines
    reflines = gpd.read_file('data/ref_lines.shp')

    fts = 10

    for j, d in enumerate(demyrs): 
        if d != 2018:   
            rf = reflines.plot(ax=ax[j], colors='k', linewidth=0.5)
            ax[j].annotate('P0', xy=(51500, 189110),  xycoords='data',
                       fontsize=fts,
                       horizontalalignment='right', verticalalignment='top',)
            ax[j].annotate('P1', xy=(51600, 188980),  xycoords='data',
                       fontsize=fts,
                       horizontalalignment='right', verticalalignment='top',)
            ax[j].annotate('P2', xy=(51720, 188760),  xycoords='data',
                       fontsize=fts,
                       horizontalalignment='right', verticalalignment='top',)
            ax[j].annotate('P3', xy=(51910, 188640),  xycoords='data',
                       fontsize=fts,
                       horizontalalignment='right', verticalalignment='top',)
        if d == 2018:
            rf01 = reflines.loc[(reflines.col1 == 'L0') | (reflines.col1 == 'L1')]
            rf01.plot(ax=ax[j], colors='k', linewidth=0.5)
            ax[j].annotate('P0', xy=(51500, 189110),  xycoords='data',
                       fontsize=fts,
                       horizontalalignment='right', verticalalignment='top',)
            ax[j].annotate('P1', xy=(51600, 188980),  xycoords='data',
                       fontsize=fts,
                       horizontalalignment='right', verticalalignment='top',)

    #fig.subplots_adjust(hspace=0, wspace=0, bottom=0.04, top=0.99, right=0.99, left=0.06)
    fig.subplots_adjust(hspace=0, wspace=0, bottom=0.04, top=0.99, right=0.99, left=0.1)
    
    # cax.set_xlabel(party)

    fig.savefig('figs/elevchange.png', dpi=300)


# overview plot
def overview(Lines, demsHEK, ref_lines):
    # load DEMs:
    el, hs, BBox = getdata(demsHEK, 9)
    # print(type(np.asarray(el)))
    ds = xr.open_rasterio(demsHEK)

    X, Y = np.meshgrid(ds[9].x.values, ds[9].y.values)
    #dem = rd.LoadGDAL('/Users/leahartl/Desktop/HEK/aligned-resam_2017.tif')
    slope = rd.TerrainAttribute(rd.rdarray(el, no_data=0), attrib='slope_degrees')

    # load centerline:
    ln = gpd.read_file('data/centerline_1.shp')  
    ln.to_crs(epsg=31254, inplace=True)

    ref_lines.to_crs(epsg=31254, inplace=True)
    
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))#, sharex=True, sharey=True)

    bounds = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45])
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

    s = ax.imshow(slope, cmap="viridis", norm=norm, alpha=0.9, extent=BBox)
    ax.imshow(hs, cmap="Greys", alpha=0.2, extent=BBox)

    CS = ax.contour(X, Y, ds[9].values, colors='k', levels=list(range(2000, 3000, 100)), linewidths=0.5)
    ax.clabel(CS, CS.levels, inline=True, fontsize=10)
            
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(s, label='slope angle', cax=cax) #Adding colorbar with label
    scalebar = ScaleBar(1)  # 1 pixel = 1 meter
    ax.add_artist(scalebar)
    msize=2
    Lines['L0']['gdf'][Lines['L0']['gdf']['Date']>'2015-01-01'].plot(ax=ax, markersize=msize, label='L0', color='k')
    Lines['L1']['gdf'][Lines['L1']['gdf']['Date']>'2015-01-01'].plot(ax=ax, markersize=msize, label='L1', color='k')
    Lines['L2']['gdf'][Lines['L2']['gdf']['Date']>'2015-01-01'].plot(ax=ax, markersize=msize, label='L2', color='k')
    Lines['L3']['gdf'][Lines['L3']['gdf']['Date']>'2015-01-01'].plot(ax=ax, markersize=msize, label='L3', color='k')
    
    # comment out this line to remove longditudinal profile
    Lines['L4']['gdf'][Lines['L4']['gdf']['Date']>'2015-01-01'].plot(ax=ax, markersize=msize, label='L4', color='k')

    fts = 10

    ln.plot(ax=ax, color='k', linewidth=0.8)
    ref_lines.plot(ax=ax, color='k', linewidth=0.8)
    ax.annotate('P0', xy=(51500, 189110),  xycoords='data',
                fontsize=fts,
                horizontalalignment='right', verticalalignment='top',)
    ax.annotate('P1', xy=(51600, 188980),  xycoords='data',
                fontsize=fts,
                horizontalalignment='right', verticalalignment='top',)
    ax.annotate('P2', xy=(51720, 188760),  xycoords='data',
                fontsize=fts,
                horizontalalignment='right', verticalalignment='top',)
    ax.annotate('P3', xy=(51910, 188640),  xycoords='data',
                fontsize=fts,
                horizontalalignment='right', verticalalignment='top',)

    ax.set_xlim([51000, 52500])
    ax.set_ylim([187980, 189320])
    # ax.legend()
    fig.savefig('figs/overview.png', dpi=150)
    # plt.show()


# compound plot of map and strain rates at P1 and P1
def figMap(Lines, bg, demsHEK):

    dfs = pd.DataFrame(columns=[2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021])
    L0 = Lines['L0']['gdf']
    # L0 = L0.loc[L0.year>2011]
    L1 = Lines['L1']['gdf']
    # L1 = L1.loc[L1.year>2011]
    L2 = Lines['L2']['gdf']
    # L2 = L2.loc[L2.year>2011]
    L3 = Lines['L3']['gdf']
    # L3 = L3.loc[L3.year>2011]

    for i, y in enumerate(dfs.columns):#[1:]):
        year1 = dfs.columns[i]
        L1_1 = L1.loc[L1.year == year1]
        L2_1 = L2.loc[L2.year == year1]

        y1 = pd.DataFrame(columns=['L1Stone', 'L2Stone', 'L1x', 'L1y', 'L1z', 'L2x', 'L2y', 'L2z'])
        y1['L1Stone'] = L1_1.Stone.values
        y1['L2Stone'] = L2_1.Stone.values[1:-1]

        y1['L1x'] = L1_1.x.values
        y1['L1y'] = L1_1.y.values
        y1['L1z'] = L1_1.z.values
        y1['L1v'] = L1_1.dxyz_a.values

        y1['L2x'] = L2_1.x.values[1:-1]
        y1['L2y'] = L2_1.y.values[1:-1]
        y1['L2z'] = L2_1.z.values[1:-1]
        y1['L2v'] = L2_1.dxyz_a.values[1:-1]

        y1['strain'] = (y1['L1v'] - y1['L2v'] ) / np.sqrt((y1['L1x']-y1['L2x'])**2 + (y1['L1y']-y1['L2y'])**2+ (y1['L1z']-y1['L2z'])**2)
        dfs[y] = y1['strain'].values
        print(y1)

    l2stns = L2.loc[L2.year == 2009].Stone.values[1:-1]

    fig, ax = plt.subplots(1, 2, figsize=(9, 6))
    ax = ax.flatten()
    el, hs, BBox = getdata(demsHEK, 13) #13--->2021 with i+1 in getdata...
    ax[0].imshow(hs, cmap="Greys", alpha=0.9, extent=BBox)
    ax[0].grid('both')

    L1.loc[L1.year>=2009].plot(ax=ax[0], markersize=1, color='red')#label='L1')
    L2.loc[L2.year>=2009].plot(ax=ax[0], markersize=1, color='k')
    L2.loc[(L2.year>=2009)&(L2.Stone>1)&(L2.Stone<12)].plot(ax=ax[0], markersize=4, color='red')

    ax[0].set_xlim([51000, 51800])
    ax[0].set_ylim([188500, 189200])
    ax[0].set_xticks([51200, 51400, 51600])
    ax[0].set_yticks([188500, 188700, 188900, 189100])

    clr = cm.plasma(np.linspace(0, 1, len(dfs.columns)+1))
    
    for i, c in enumerate(dfs.columns):
        # alternate marker types so that years are easier to visually distinguish
        if i % 2 == 0:
            ax[1].scatter(dfs.index+1, dfs[c].values, label=c, color=clr[i], marker='o')
        else:
            ax[1].scatter(dfs.index+1, dfs[c].values, label=c, color=clr[i], marker='^')
        ax[1].legend()
    ax[1].set_ylabel('strain rate (a$^{−1}$)')
    ax[1].set_xlabel('block nr.')
    ax[1].legend()
    ax[1].grid('both')

    bbox = dict(boxstyle="round", fc="0.99")
    ax[0].annotate('a)', xy=(0.95, 0.9),  xycoords='axes fraction',
                   fontsize=10, horizontalalignment='right', verticalalignment='top', bbox=bbox,)
    ax[1].annotate('b)', xy=(0.95, 0.9),  xycoords='axes fraction',
                   fontsize=10, horizontalalignment='right', verticalalignment='top', bbox=bbox,)

    plt.tight_layout()
    fig.savefig('figs/figMapStrain2.png', dpi=300)


def plotRotation(demsHEK, dhHEK):
    # load shape file with angles (rot_roll angle in file is the picth angle, which we want... misnamed in the file.)
    data = gpd.read_file('/Users/leahartl/Desktop/HEK/magnus2/ICP_Match_Results_shifted.shp')
    data.to_crs(epsg=31254, inplace=True)

    # separate the years.
    d2019 = data.loc[data['ID']==2019]
    d2020 = data.loc[data['ID']==2020]
    d2021 = data.loc[data['ID']==2021]

    # start plot
    al = 0.9
    fig, ax = plt.subplots(1, 3, figsize=(10, 6), sharex=True, sharey=True)
    ax = ax.flatten()
    divider = make_axes_locatable(ax[-2])
    # set DEMs to be plotted
    demyrs = [1971, 1977, 1990, 1997, 2006, 2009, 2010, 2011, 2017, 2018, 2019, 2020, 2021]
    ttlyrs = ['1953-71', '1971-77', '1977-90', '1990-97', '1997-06', '2006-09', '2009-10', '2010-11',
              '2011-17', '2017-18', '2018-19', '2019-20', '2020-21']
    # Open the DEM with Rasterio, make and plot hillshades with earthpy
    for i, d in enumerate([11, 12, 13]):
        print(d)
        el, hs, BBox = getdata(demsHEK, d)
        im = ax[i].imshow(hs, cmap="Greys", alpha=al, extent=BBox)
        ax[i].grid('both')

    #plot elevation change rasters (careful with numbering, DDSM stack has one less layer than the DSM stack)    
    for i, d in enumerate([10, 11, 12]):
        # print(d)
        el, BBox2 = getdata2(dhHEK, d)
        alphas = np.ones(el.shape)*0.6
        alphas[(el < 0.2) & (el>-0.2)] = 0.2
        alphas[(el < -100)] = 0

        im2 = ax[i].imshow(el, cmap='RdBu', alpha=alphas, extent=BBox2, vmin=-4, vmax=+4)
        ax[i].set_title(ttlyrs[d-1])
        ax[i].grid('both')

    cmrot = cm.get_cmap('Greys', 5)
    s = 24
    # cmrot = 'PiYG'
    
    d2020 = d2020.reset_index()
    d2021 = d2021.reset_index()

    d=d2019.plot(column='ROT_ROLL_me', cmap=cmrot, vmin=-5, vmax=5, ax=ax[0], s=s, edgecolors='k')
    d2020.plot(column='ROT_ROLL_me', cmap=cmrot, vmin=-5, vmax=5, ax=ax[1], s=s, edgecolors='k')
    d2021.plot(column='ROT_ROLL_me', cmap=cmrot, vmin=-5, vmax=5, ax=ax[2], s=s, edgecolors='k')

    fig.subplots_adjust(right=0.9)
    cbar_ax1 = fig.add_axes([0.91, 0.3, 0.02, 0.4])
    cbim = fig.colorbar(cm.ScalarMappable(cmap='RdBu', norm=colors.Normalize(vmin=-4, vmax=4)), cax=cbar_ax1, orientation='vertical')
    cbim.set_label('elevation change (m)')

    ax[0].set_xlim([51150, 51600])
    ax[0].set_ylim([188520, 189150])
    ax[0].set_xticks([51300, 51500])


    # scalebar = ScaleBar(1) # 1 pixel = 1 meter
    cbar_ax = fig.add_axes([0.5, 0.1, 0.3, 0.02])

    cbpoints = fig.colorbar(cm.ScalarMappable(cmap=cmrot, norm=colors.Normalize(vmin=-5, vmax=5)), cax=cbar_ax, orientation='horizontal')
    cbpoints.set_label('pitch angle (°)')

    plt.subplots_adjust(hspace=0.01, wspace=0, bottom=0, top=1)

    fig.savefig('figs/elevchangeRotation3.png', dpi=300)



