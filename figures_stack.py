#! /usr/bin/env python3
import numpy as np
import pandas as pd
# import fiona
import geopandas as gpd
import rasterio
import matplotlib.pyplot as plt
from matplotlib import cm
# from scipy.interpolate import UnivariateSpline
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset, InsetPosition
# import xarray as xr
# import rioxarray
from scipy import signal
from matplotlib.gridspec import GridSpec


# ---- some helper functions -----#
def getdata(f, i):
    with rasterio.open(f) as src:
        el = src.read(i-1)
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


def BCF(co, alpha, H2):
    n = 2.1
    rho = 1500 #2000
    o = 25 # alpha
    tau = 10*1000
    y = 0.06
    g = 9.81
    H1 = 100  *1000/(rho*g*np.sin(alpha*np.pi / 180.))
    # print('H', H)
    F1 = (n+1)/y
    F2 = tau + rho * g * H2 * np.cos(alpha* np.pi / 180.)*np.tan(o*np.pi / 180.)
    F2_c = tau + rho * g * H1 * np.cos(alpha* np.pi / 180.)*np.tan(o*np.pi / 180.)

    F3 = (rho*g*np.sin(alpha * np.pi / 180.))

    F4 = H2**-(n+1)
    F4_c = H1**-(n+1)

    BCF = co * F1 * (F2/F3)**n *F4
    BCF_c = co * F1 * (F2_c/F3)**n *F4_c
    return(BCF, H1)


def extractFlowline(fname, xs, ys):
    with rasterio.open(fname) as src:
        out = np.asarray([smpl[0] for smpl in src.sample(zip(xs, ys))])
        out[out < 0] = np.nan
    return (out)


def extractFlowlineStack(fname, xs, ys, yrs):
    # demyrs = [1953, 1971, 1977, 1990, 1997, 2006, 2009, 2010, 2017, 2018, 2019, 2020, 2021]

    out = pd.DataFrame(columns=yrs)
    with rasterio.open(fname) as src:
        for i, c in enumerate(out):
            out[c] = np.asarray([smpl[i] for smpl in src.sample(zip(xs, ys))])
            #out[c][out[c] < 0] = np.nan
    return (out)


def collectFlowline(demsHEK):
    demyrs = [1953, 1971, 1977, 1990, 1997, 2006, 2009, 2010, 2011, 2017, 2018, 2019, 2020, 2021]
    # exctract points from flowline
    ln = gpd.read_file('data/centerline_1.shp') 
    ln.to_crs(epsg=31254, inplace=True) 
    line =ln.geometry[0]
    distance_delta = 1
    distances = np.arange(0, line.length, distance_delta)
    points = [line.interpolate(distance) for distance in distances] + [line.boundary[1]]
    xs = [point.x for point in points]
    ys = [point.y for point in points]
    d = pd.DataFrame(columns=['x', 'y'])
    d['x'] = xs
    d['y'] = ys
    d['dxy'] = np.sqrt((d['x'].diff()**2 + d['y'].diff()**2))
    # extract data from DEMs for each point along flowline
    flowdata = pd.DataFrame(index=d['dxy'].values)
    
    #for i, d in enumerate(demsHEK):
    flowdata = extractFlowlineStack(demsHEK, xs, ys, demyrs)
    # print(flowdata)
    # stop
    ug = extractFlowline('data/DEMs/test_aligned-untergrund1_31254_GEO.tif', xs, ys)

    # hardcoded removal of artefacts in dems.
    ug[1600:] = np.nan

    
    ug[0:266] = flowdata[1953].iloc[0:266]
    # # ug[0:152] = flowdata[1971].iloc[0:152]
    ug[0:200] = flowdata[1971].iloc[0:200]
    flowdata.where(flowdata > 0, np.nan, inplace=True)

    return(flowdata, ug)


def getcenterline():
    ln = gpd.read_file('data/centerline_1.shp')
    ln.to_crs(epsg=31254, inplace=True)  
    line =ln.geometry[0]
    distance_delta = 1
    distances = np.arange(0, line.length, distance_delta)
    points = [line.interpolate(distance) for distance in distances] + [line.boundary[1]]
    xs = [point.x for point in points]
    ys = [point.y for point in points]
    d = pd.DataFrame(columns=['x', 'y'])
    d['x'] = xs
    d['y'] = ys
    d['dxy'] = np.sqrt((d['x'].diff()**2 + d['y'].diff()**2))
    # print(d['dxy'])
    return (d, xs, ys)


def getDepthAngle(demsHEK, demyrs):
    # demyrs = [1953, 1971, 1977, 1990, 1997, 2006, 2009, 2010, 2011, 2017, 2018, 2019, 2020, 2021]
    # get points along line
    d, xs, ys = getcenterline()

    # get data for the points
    demdata, ug = collectFlowline(demsHEK)
    # remove negative no data values
    demdata.where(demdata > 0, np.nan, inplace=True)

    # generate dataframe of depth along flowline (subtract bedrock, smooth it a bit)
    depth = demdata.sub(pd.Series(ug), axis=0)
    depth = depth.rolling(30).mean()

    demdata.index = d['dxy']
    demdata.index.values[0] = 0

    depth.index = d['dxy']
    depth.index.values[0] = 0

    # generate dataframe of angles along flowline
    a = demdata.diff().div(demdata.index.values, axis=0)
    a = a.apply(np.arctan) * (180/np.pi)
    a2 = a.rolling(50).mean()
    a2.where(a2 > 0, np.nan, inplace=True)

    return (demdata, ug, depth, a2)


def BCF_TS(meanV, ds, demsHEK, demyrs):
    #demyrs = [1953, 1971, 1977, 1990, 1997, 2006, 2009, 2010, 2011, 2017, 2018, 2019, 2020, 2021]
    demdata, ug, depth, a2 = getDepthAngle(demsHEK, demyrs)

    dat2 = pd.read_csv('data/mean_velocities_block_monitoring.txt', delimiter='\t')
    dat2.index = pd.to_datetime(dat2['date1'], format="%d.%m.%Y").dt.year
    dat2.index.name = 'year'

    cs = a2.index.values.cumsum()
    a2.index = cs
    depth.index = cs

    ds = [round(elem, 0) for elem in ds] 
    # extract values for profile positions
    val = ds[0]
    s = pd.Series(a2.index.values).sub(val).abs().idxmin()

    p_d = pd.DataFrame(columns=['d0', 'd1', 'd2', 'd3'], index=demyrs)
    for i, p in enumerate(p_d.columns):
        p_d[p] = depth.loc[ds[i]:ds[i]+1, :].values.flatten().tolist()

    p_d = p_d.reindex(list(range(p_d.index.min(), p_d.index.max()+1)), fill_value=np.nan)
    
    meanV.set_index('year', drop=True, inplace=True)

    # hardcoded angles, produced in "proc_rasters.py"
    a = pd.read_csv('angles.csv').set_index('year')

    temp = pd.concat([meanV, a], axis=1, join="outer")

    result1 = pd.concat([temp, p_d], axis=1, join="inner")
    result2 = pd.concat([result1, dat2], axis=1, join="outer")
    result2.sort_index(inplace=True)
    result1.sort_index(inplace=True)

    result2.drop(columns=['P0', 'P1', 'P2', 'P3'], inplace=True)
    result2.dropna(how='all', inplace=True)
    # print(result2.head(30))

    result = result1.interpolate(method='linear', limit_direction='forward', axis=0)
    # print(result.head())
    fig, ax = plt.subplots(4, 1, figsize=(9, 7), sharex=True, sharey=True)
    ax = ax.flatten()

    Bdf = pd.DataFrame(index=result.index, columns=[0, 1, 2, 3])
    Bdf_DSM = pd.DataFrame(index=result2.index, columns=[0, 1, 2, 3])

    for p in [0, 1, 2, 3]:
        B, h = BCF(result['P'+str(p)].values, result['a'+str(p)].values, result['d'+str(p)].values)
        Bdf[p] = B

        B2, h2 = BCF(result2['vel_p'+str(p)].values, result2['a'+str(p)].values, result2['d'+str(p)].values)
        Bdf_DSM[p] = B2

        ax[p].plot(result.index.values, B, label='BCF', linestyle='-', color='k')

        ax[p].set_ylabel('BCF')
        ax[p].grid('both')
        ax[p].set_title('P'+str(p))
        ax0 = ax[p].twinx()
        color = 'grey'
        ax0.set_ylabel('mean vel. (m/a)', color=color)  
        ax0.tick_params(axis='y', labelcolor=color)
        ax0.plot(result.index.values, result['P'+str(p)], label='mean vel. (m/a)', linestyle='-', color=color)
        ax0.set_ylim([0, 14])


    # ask matplotlib for the plotted objects and their labels
    lines, labels = ax[0].get_legend_handles_labels()
    lines2, labels2 = ax0.get_legend_handles_labels()
    ax[0].legend(lines + lines2, labels + labels2, loc=0)
    # fig.savefig('figs/BCF_subpots.png')

    # fig2, ax2 = plt.subplots(1, 1, figsize=(10, 6))
    # # ax2 = ax2.flatten()

    # for p in [0, 1, 2, 3]:
    #     # B, h = BCF(result['P'+str(p)].values, aa[p], hh[p])
    #     B, h = BCF(result['P'+str(p)].values, result['a'+str(p)].values, result['d'+str(p)].values)
    #     ax2.plot(result.index.values, B, label='BCF P'+str(p), linestyle='-')#, color='k')
    # ax2.set_ylim([0, 50])

    #     # ax[p].axhline(y=1, color='k', linestyle='--')
    # ax2.set_ylabel('BCF')
    # ax2.grid('both')
    # # ax2.set_title('P'+str(p))
    # ax0 = ax2.twinx()
    # color = 'grey'
    # ax0.set_ylabel('mean vel. (m/a)', color=color)  
    # ax0.tick_params(axis='y', labelcolor=color)
    # for p in [0, 1, 2, 3]:
    #     ax0.plot(result.index.values, result['P'+str(p)], label='mean vel. P'+str(p), linestyle='--', linewidth=1)#, color=color)
    # ax0.set_ylim([0, 14])


    # # ask matplotlib for the plotted objects and their labels
    # lines, labels = ax2.get_legend_handles_labels()
    # lines2, labels2 = ax0.get_legend_handles_labels()
    # ax2.legend(lines + lines2, labels + labels2, loc=0)
    # fig2.savefig('figs/BCF_profiles.png')
    return(Bdf, Bdf_DSM)


# --- plotting functions for paper figures start here ----# 
# Velocity along flowline, data extracted from disp rasters (supplement)
def FlowlineVel(demsHEK, disp, ds, dispyrs):
    # exctract points from flowline
    d, xs, ys = getcenterline()

    # extract data from disp rasters for each point along flowline
    flowdataD = extractFlowlineStack(disp, xs, ys, dispyrs)

    flowdataD.index = d['dxy'].cumsum()
    flowdataD.index.values[0] = 0
    # remove negative outliers if any
    flowdataD[flowdataD < -100] = np.nan

    # start plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 7), sharex=True)

    interesting_keys = [1971, 1977, 1990, 1997, 2006, 2009, 2010, 2011, 2017, 2018, 2019, 2020, 2021]
    lbls = ['1953-71', '1971-77', '1977-90', '1990-97', '1997-06', '2006-09', '2009-10', '2010-11',
            '2011-17', '2017-18', '2018-19', '2019-20', '2020-21']
    # interesting_keys = [2010, 2011, 2017, 2018, 2019, 2020, 2021]
    clr = cm.plasma(np.linspace(0, 1, len(interesting_keys)+1))

    #remove artefacts at the edges of the ULS data
    flowdataD.loc[900:, 2021] = np.nan
    flowdataD.loc[650:, 2019] = np.nan
    flowdataD.loc[690:, 2018] = np.nan

    for i, c in enumerate(interesting_keys):
        ax.plot(flowdataD.index.values, flowdataD[c].values, label=lbls[i], color=clr[i], linewidth=1)

    ax.grid('both')
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='center right')
    ax.set_ylabel('velocity (m/a)')
    ax.set_xlabel('distance along flow line (m)')
    ax.set_xlim([100, 1000])
    ax.set_ylim([0, 30])

    plt.tight_layout()
    fig.savefig('figs/FlowlineVel.png')


def SingleBlocks(Lines, pts):
    fig, axs = plt.subplots(3, 2, figsize=(10, 7), sharex=False, sharey=False,
                            gridspec_kw={'width_ratios': [1, 1], 'height_ratios': [1, 1, 1]})
    #ax = ax.flatten()
    gs = axs[1, 1].get_gridspec()
    for ax in axs[-1, :]:
        ax.remove()
    axbig = fig.add_subplot(gs[-1, :])

    ax = axs.flatten()
    clrs = cm.get_cmap('viridis', 6).colors

    Lines['L0']['dat']['St'] = Lines['L0']['dat']['Stone']
    Lines['L0']['dat'] = Lines['L0']['dat'].dropna()

    Lines['L1']['dat']['St'] = Lines['L1']['dat']['Stone'].str.replace(r'\D', '')
    Lines['L1']['dat'] = Lines['L1']['dat'].dropna()

    Lines['L2']['dat']['St'] = Lines['L2']['dat']['Stone']
    Lines['L2']['dat'] = Lines['L2']['dat'].dropna()

    Lines['L3']['dat']['St'] = Lines['L3']['dat']['Stone']
    Lines['L3']['dat'] = Lines['L3']['dat'].dropna()

    L4 = Lines['L4']['gdf']
    # print(L4)

    yr = Lines['L1']['dat']['Date'].dt.year.unique()

    interesting_keys = ('L0', 'L1', 'L2', 'L3')
    subdict = {x: Lines[x] for x in interesting_keys if x in Lines}

    for i, L in enumerate(subdict):
        # print(i, L)
        for j, y in enumerate(yr[yr > 2015]):
            temp = Lines[L]['dat'].loc[Lines[L]['dat']['Date'].dt.year == y]
            pts2 = pts[pts.p == 'P'+str(i)+'S']
            ## something is messed up here with x y column labeling; xy are reversed in either in pts or Lines gdf
            dis = np.sqrt(((pts2['x'].values-temp['y'])**2 + (pts2['y'].values-temp['x'])**2))
            ax[i].scatter(dis, temp['dxyz_a'], label=str(y), color=clrs[j])#, cmap='viridis')
        ax[i].grid(axis='both')
        ax[i].set_title('P'+str(i))
        ax[i].set_ylabel('velocity (m/a)')
        ax[i].set_xlabel('distance from reference point (m)')

    for s in [5.0, 6.0, 7.0, 8.0, 9.0, 10.0]: # only use stones available in all years 2016-2021
        for j, y in enumerate(yr[yr > 2015]):
            temp = Lines['L4']['dat'].loc[Lines['L4']['dat']['Date'].dt.year == y]
            axbig.scatter(temp['z'], temp['dxyz_a'], color=clrs[j])
            
    ax[0].legend()
    ax[0].set_ylim([0, 22])
    ax[1].set_ylim([0, 22])
    ax[2].set_ylim([0, 22])
    ax[3].set_ylim([0, 22])
    axbig.set_ylim([0, 22])
    axbig.grid('both')
    axbig.set_ylabel('velocity (m/a)')
    axbig.set_xlabel('elevation (m.a.s.l.)')
    axbig.set_title('P long.')
    plt.tight_layout()
    fig.savefig('figs/SingleBlocks.png', dpi=300)


def timeseriesBoth4(meanV, Bdf, Bdf_DSM, demyrs):
    # 3 subplots, time series of mean block profile vel, DSM vel, BCF; as horizontal lines.

    # load data and prepare for plotting horiz. lines.
    dat2 = pd.read_csv('data/mean_velocities_block_monitoring.txt', delimiter='\t')
    dat2.index = pd.to_datetime(dat2['date1'], format="%d.%m.%Y").dt.year+1
    dat2.index.name = 'year'
    dat3 = dat2.reindex(np.arange(1953, 2022), method='ffill')
    dat3.index.name = 'year'
    dat3.drop(columns=['date1'], inplace=True)
    
    dat2['year1'] = pd.to_datetime(dat2['date1'], format="%d.%m.%Y").dt.year
    dat2['year2'] = pd.to_datetime(dat2['date1'].shift(-1), format="%d.%m.%Y").dt.year

    meanV['year2'] = meanV.index.values
    meanV['year1'] = meanV['year2'].shift(1)
    meanV.loc[1952, 'year1'] = 1951

    Bdf['year2'] = Bdf.index.values
    Bdf['year1'] = Bdf['year2'].shift(1)

    Bdf_DSM['year2'] = Bdf_DSM.index.values
    Bdf_DSM['year1'] = Bdf_DSM['year2'].shift(-1)
    Bdf_DSM.drop([2021], inplace=True)

    # start figure
    fig, ax = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
    ax = ax.flatten()

    # set colors
    # cl = ['#cc79a7', '#0072b2', '#f0e442', '#009e73']
    cl = ['#BB5566', '#004488', '#228833', '#DDAA33']
    for i, p in enumerate([0, 1, 2, 3]):
        ax[0].hlines(y=meanV['P'+str(p)], xmin=meanV.year1, xmax=meanV.year2, color=cl[i], label='P' + str(p))

    for i, p in enumerate(dat3.columns):
        ax[1].hlines(y=dat2[p], xmin=dat2.year1, xmax=dat2.year2, color=cl[i])

    ax[1].set_ylabel('mean velocity (m/a)')
    ax[0].grid('both')
    ax[1].grid('both')
    ax[0].legend(loc='upper left')
    ax[0].set_ylabel('mean velocity (m/a)')
    ax[0].set_title('Block profiles')
    ax[1].set_title('DSM-based image correlation')

    for y in demyrs:
        ax[1].grid('both')
        ax[1].vlines(x=y, ymin=0, ymax=ax[1].get_ylim()[1], colors='grey', linewidth=1, linestyle='--',)
    
    ax[0].set_ylim([0, 17])
    ax[1].set_ylim([0, 17])

    for i, p in enumerate([0, 1, 2, 3]):
        ax[2].hlines(y=Bdf[p], xmin=Bdf.year1, xmax=Bdf.year2, color=cl[i])
        # careful: years shifted the other way in BDF DSM dataframe
        ax[2].hlines(y=Bdf_DSM[p], xmin=Bdf_DSM.year2, xmax=Bdf_DSM.year1, color=cl[i], linestyle='--', linewidth=1)
    ax[2].set_ylim([0, 44])
    ax[2].set_xlim([1952, 2022])
    ax[2].set_ylabel('BCF')
    ax[2].grid('both')
    ax[2].set_title('Bulk creep factor')

    bbox = dict(boxstyle="round", fc="0.99")
    a0 = ax[0].annotate('a)', xy=(0.15, 0.9),  xycoords='axes fraction',
                        fontsize=10, horizontalalignment='right', verticalalignment='top',bbox=bbox,)
    a1 = ax[1].annotate('b)', xy=(0.05, 0.9),  xycoords='axes fraction',
                        fontsize=10, horizontalalignment='right', verticalalignment='top',bbox=bbox,)
    a2 = ax[2].annotate('c)', xy=(0.05, 0.9),  xycoords='axes fraction',
                        fontsize=10, horizontalalignment='right', verticalalignment='top',bbox=bbox,)
    plt.tight_layout()
    fig.savefig('figs/TS_BCF_meanV_4.png', dpi=300)


def FlowlineAll_inset(demsHEK, ds):
    ## main flowline plot with inset axies
    # exctract points from flowline
    d, xs, ys = getcenterline()

    flowdata, ug = collectFlowline(demsHEK)
    ug[0:266] = np.nan
    flowdata.index = d['dxy']
    flowdata.index.values[0] = 0

    # start figure
    fig, ax = plt.subplots(figsize=(12, 7))
    # plot bedrock
    ax.plot(flowdata.index.values.cumsum(), ug, label='Bedrock', color='k', linestyle='--', linewidth=0.5)

    interesting_keys = [1953, 1971, 1977, 1990, 1997, 2006, 2009, 2010, 2011, 2017, 2018, 2019, 2020, 2021]
    # set colormap
    clr = cm.magma(np.linspace(0, 1, len(interesting_keys)+1))
    flowdata.where(flowdata > 0, np.nan, inplace=True)

    # plot lines
    for i, c in enumerate(flowdata[interesting_keys].columns):
        if i % 2 == 0:
            ax.plot(flowdata.index.values.cumsum(), flowdata[c].values, label=c, color=clr[i], linewidth=0.7, linestyle='--')
        else:
            ax.plot(flowdata.index.values.cumsum(), flowdata[c].values, label=c, color=clr[i], linewidth=0.7)

    # add vertical lines and annotations
    elv = [2460, 2555, 2655, 2710]
    ymn = [2400, 2500, 2600, 2660]
    ymx = [2460, 2550, 2650, 2710]
    for i, d in enumerate(ds):
        ax.vlines(x=ds[i], ymin=ymn[i], ymax=ymx[i], colors='k', linewidth=2)
        ax.annotate('P'+str(i), xy=(ds[i], elv[i]),  xycoords='data',
                    fontsize=10, horizontalalignment='right', verticalalignment='bottom',)

    ax.annotate('A', xy=(530, 2568),  xycoords='data',
                fontsize=16, weight='bold', horizontalalignment='right', verticalalignment='top',)
    ax.annotate('B', xy=(386, 2498),  xycoords='data',
                fontsize=16, weight='bold', horizontalalignment='right', verticalalignment='top',)

    ax.grid('both')
    ax.legend()
    ax.set_ylabel('elevation (m.a.s.l.)')
    ax.set_xlabel('distance along central flowline (m)')
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')
    ax.set_ylim([2300, 2750])
    ax.set_xlim([0, 1000])

    # set up inset axes
    axins = zoomed_inset_axes(ax, zoom=2.5, loc='lower right')
    axins.set_xlim(360, 570)
    axins.set_ylim(2500, 2590)
    axins.grid('both')
    axins.plot(flowdata.index.values.cumsum(), ug, label='Bedrock', color='k', linestyle='--', linewidth=0.5)
    for i, c in enumerate(flowdata[interesting_keys].columns):
        if i % 2 == 0:
            axins.plot(flowdata.index.values.cumsum(), flowdata[c].values, label=c, color=clr[i], linewidth=0.7, linestyle='--')
        else:
            axins.plot(flowdata.index.values.cumsum(), flowdata[c].values, label=c, color=clr[i], linewidth=0.7)

    # draw a bbox of the region of the inset axes in the parent axes and
    # connecting lines between the bbox and the inset axes area
    mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.5")
    # right, bottom, width, height
    ip = InsetPosition(ax, [.5, -0.1, 0.54, 0.52])
    axins.set_axes_locator(ip)
    axins.annotate('1st destab. sign', xy=(515, 2573),  xycoords='data', xytext=(528, 2568),
                   fontsize=12, horizontalalignment='left', verticalalignment='top',
                   arrowprops=dict(facecolor='black', shrink=0.05),)
    axins.annotate('major scarp 1st phase', xy=(400, 2509),  xycoords='data', xytext=(490, 2505),
                   fontsize=12, horizontalalignment='right', verticalalignment='top',
                   arrowprops=dict(facecolor='black', shrink=0.05),)
    axins.annotate('major scarp 2nd phase', xy=(507, 2562),  xycoords='data', xytext=(520, 2545),
                   fontsize=12, horizontalalignment='center', verticalalignment='top',
                   arrowprops=dict(facecolor='black', shrink=0.05),)

    axins.yaxis.get_major_locator().set_params(nbins=4)
    axins.xaxis.get_major_locator().set_params(nbins=4)

    fig.savefig('figs/flowline_all_inset.png', dpi=300)


def FlowlineScarp(demsHEK, ds):
    # exctract points from flowline
    d, xs, ys = getcenterline()

    flowdata, ug = collectFlowline(demsHEK)
    ug[0:266] = np.nan

    flowdata.index = d['dxy'].cumsum()
    flowdata.index.values[0] = 0

    # define area around scarp of interest
    flowdata = flowdata.loc[flowdata.index[285:420]]

    # make copy of data an detrend it
    flowdata2 = flowdata.copy()
    for c in flowdata.columns:
        flowdata2[c] = signal.detrend(flowdata[c])

    # remove issue in 1977 line
    flowdata2.loc[0:340, 1977] = np.nan

    # get minimum
    ixmn = flowdata2.idxmin()
    mn = flowdata2.min()

    # make another copy, set value above minimum indeces to nan to get first max downslope of minimum in detrended data
    flowdata3 = flowdata2.copy()
    for c in flowdata3:
        flowdata3[c] = flowdata3[c].loc[flowdata3.index < flowdata2[c].idxmin()]

    ixmx = flowdata3.idxmax()
    mx = flowdata3.max()

    # same but for original, not detrended daata
    flowdata3_org = flowdata.copy()
    for c in flowdata3_org:
        flowdata3_org[c] = flowdata3_org[c].loc[flowdata3_org.index < flowdata2[c].idxmin()]

    flowdata4 = flowdata.copy()
    last = pd.DataFrame(columns=flowdata3.columns, index=[1])
    val = pd.DataFrame(columns=flowdata3.columns, index=[1])

    # use shifts of original data to pick out maxima
    for c in flowdata4:
        flowdata4[c] = flowdata3_org[c][(flowdata3_org[c].shift(2) < flowdata3_org[c]) & (flowdata3_org[c].shift(-1) < flowdata3_org[c])]
        last[c] = flowdata4[c].last_valid_index()

    # arrange it all into a dataframe for  plotting
    flowdata4n = flowdata4.dropna(how='all')
    for c in flowdata4:
        tmp = flowdata4[c].dropna()
        if len(tmp) > 0:
            val[c] = tmp.values[-1]
        else:
            val[c] = np.nan


    # manually fix 2021 values
    last[2021] = ixmx[2021]
    val[2021] = mx[2021]

    df = pd.DataFrame(index=mn.index)
    df['distmx'] = last.loc[1, :]
    df['distmn'] = ixmn
    df['mx'] = val.loc[1, :]
    df['mn'] = mn
    df.loc[1953, 'mn'] = np.nan
    df['magn'] = df['mx'] - df['mn']
    df['length'] = df['distmn'] - df['distmx']

    fll = pd.DataFrame(index=df.index.values[1:])
    for n in df['distmx'].index.values[1:]:
        fll.loc[n, 'max'] = flowdata.loc[df.loc[n, 'distmx'], n]
        fll.loc[n, 'mn'] = flowdata.loc[df.loc[n, 'distmn'], n]

    # manually adjust max value for 2021 round 2, the shift is not catching it in this case.
    df.loc[2021, 'distmx'] = 298
    fll.loc[2021, 'max'] = 2477.22

    # start figure
    fig = plt.figure(figsize=(14, 6))
    gs = GridSpec(1, 2, left=0.20, wspace=0.14)#, width_ratios=[1, 1], height_ratios=[1, 1])
    ax1 = fig.add_subplot(gs[0, 0])
    ax4 = fig.add_subplot(gs[0, 1])

    interesting_keys = [1953, 1971, 1977, 1990, 1997, 2006, 2009, 2010, 2011, 2017, 2018, 2019, 2020, 2021]

    clr = cm.magma(np.linspace(0, 1, len(interesting_keys)+1))

    # plot lines
    for i, c in enumerate(interesting_keys):
        # alternate line styles
        if i % 2 == 0:
            ax1.plot(flowdata.index.values, flowdata[c].values, label=c, color=clr[i], linewidth=0.7, linestyle='--')
        else:
            ax1.plot(flowdata.index.values, flowdata[c].values, label=c, color=clr[i], linewidth=0.7)

    ax1.scatter(df['distmx'].values[1:], fll['max'].values, color='grey', label='high point')
    ax1.scatter(df['distmn'].values[1:], fll['mn'].values, color='red', label='low point')

    ax4.scatter(df.index, df['length'].values, label='length', color='grey')
    ax4.scatter(fll.index, fll['max'] - fll['mn'], color='k', label='height')
    
    ax1.grid('both')
    ax4.grid('both')
    ax1.set_title('surface elevation along flow line')
    ax1.legend(bbox_to_anchor=(-0.15, 0.8))
    ax1.set_ylabel('elevation (m.a.s.l.)')
    ax1.set_xlabel('distance along central flowline (m)')
    ax4.set_title('size of scarp')
    ax4.set_xlabel('years of DSM acquisition')
    ax4.set_ylabel('size (m)')
    ax4.legend()

    bbox = dict(boxstyle="round", fc="0.99")
    ax1.annotate('a)', xy=(292, 2530),  xycoords='data',
                 fontsize=10, horizontalalignment='right', verticalalignment='bottom', bbox=bbox)
    ax4.annotate('b)', xy=(1972, 12.8),  xycoords='data',
                 fontsize=10, horizontalalignment='right', verticalalignment='bottom', bbox=bbox)

    fig.savefig('figs/flowline_all_scarp.png', dpi=300)


def FlowlineTerminus(demsHEK):
    # exctract points from flowline
    d, xs, ys = getcenterline()
    
    flowdata, ug = collectFlowline(demsHEK)
    ug[0:266] = np.nan
    flowdata.index = d['dxy']
    flowdata.index.values[0] = 0

    demyrs = [1953, 1971, 1977, 1990, 1997, 2006, 2009, 2010, 2011, 2017, 2018, 2019, 2020, 2021]
    # plot only selected years so changes are easier to see.
    interesting_keys = [1953, 1971, 1977, 1990, 1997, 2006, 2017, 2021]#, 2009, 2010, 2011, 2017, 2018, 2019, 2020, 2021]
    # clr = cm.cividis(np.linspace(0, 1, len(demyrs)))
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))
    for i, c in enumerate(interesting_keys):
        clr = cm.cividis(np.linspace(0, 1, len(interesting_keys)))
        if i % 2 == 0:
            ax.plot(flowdata.index.values.cumsum(), flowdata[c], label=c, linewidth=1, color=clr[i])
        else: 
            ax.plot(flowdata.index.values.cumsum(), flowdata[c], label=c, linewidth=1, color=clr[i], linestyle='--')
    ax.plot(flowdata.index.values.cumsum(), ug, label='BR', color='k', linestyle='--', linewidth=0.5)


    ax.set_xlim([100, 570])
    ax.set_ylim([2350, 2590])

    ax.set_ylabel('elevation (m.a.s.l.)')
    ax.set_xlabel('distance along central flowline (m)')
    ax.legend()
    ax.grid('both')
    fig.savefig('figs/flowline_terminus.png', dpi=300)

