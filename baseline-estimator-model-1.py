#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#------------------------------------------------------------------------------
# PROGRAM: baseline-estimator-model-1.py
#------------------------------------------------------------------------------
# Verion 0.3
# 7 July, 2021
# Michael Taylor
# https://patternizer.github.io
# patternizer AT gmail DOT com
# michael DOT a DOT taylor AT uea DOT ac DOT uk
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# IMPORT PYTHON LIBRARIES
#------------------------------------------------------------------------------
# Dataframe libraries:
import numpy as np
import numpy.ma as ma
import itertools
import pandas as pd
import xarray as xr
import pickle
from datetime import datetime
import nc_time_axis
import cftime
import re

# Plotting libraries:
import matplotlib
import matplotlib.pyplot as plt; plt.close('all')
import matplotlib.cm as cm
from matplotlib import colors as mcol
from matplotlib.cm import ScalarMappable
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
from matplotlib.ticker import FuncFormatter
from matplotlib.collections import PolyCollection
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
import cmocean
import seaborn as sns; sns.set()

# OS libraries:
import os
import os.path
from pathlib import Path
import sys
import subprocess
from subprocess import Popen
import time

# Stats libraries:
import random
from sklearn import datasets, linear_model
from sklearn.linear_model import LinearRegression, TheilSenRegressor
from sklearn.linear_model import RANSACRegressor
from sklearn.metrics import mean_squared_error, r2_score
from statsmodels.tsa.stattools import adfuller
import statsmodels.api as sm
from statsmodels.nonparametric.smoothers_lowess import lowess

# Maths libraries:
from math import radians, cos, sin, asin, sqrt

# Datetime libraries
import cftime
import calendar 
# print(calendar.calendar(2020))
# print(calendar.month(2020,2))
# calendar.isleap(2020)
from datetime import date, time, datetime, timedelta
#today = datetime.now()
#tomorrow = today + pd.to_timedelta(1,unit='D')
#tomorrow = today + timedelta(days=1)
#birthday = datetime(1970,11,1,0,0,0).strftime('%Y-%m-%d %H:%M')
#print('Week:',today.isocalendar()[1])

# Silence library version notifications
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", message="SettingWithCopyWarning")
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# SETTINGS: 
#------------------------------------------------------------------------------

fontsize = 16
nsmooth = 60

load_glosat = True
plot_monthly = True
plot_fit = True

use_fahrenheit = True
if use_fahrenheit: 
    temperature_unit = 'F'
else:
    temperature_unit = 'C'

segment_start = pd.to_datetime('1851-01-01')
segment_end = pd.to_datetime('1900-12-01')
normal_start = pd.to_datetime('1961-01-01')
normal_end = pd.to_datetime('1990-12-01')

test_station = '744920'     # BHO
test_radius = 312           # km
n_baseline_years = 15       # Minimum number of years in baseline for calculation of normals 
n_overlap_years = 25        # Minimum number of years overlap in the segment epoch

#------------------------------------------------------------------------------
# METHODS: 
#------------------------------------------------------------------------------

def fahrenheit_to_centigrade(x):
    y = (5.0/9.0) * (x - 32.0)
    return y

def centigrade_to_fahrenheit(x):
    y = (x * (9.0/5.0)) + 32.0
    return y

def haversine(lat1, lon1, lat2, lon2):

    # CALCULATE: Haversine distance in km

    lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2]) # decimal degrees --> radians    
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    km = 6371* c # Earth radius = 6371 km
    return km

def find_nearest_lasso(lat, lon, df, radius):
    distances = df.apply(lambda row: haversine(lat, lon, row['lat'], row['lon']), axis=1)
    lasso = distances < radius
    return df.loc[lasso,:]
    
def calculate_normals_and_SEs_uncorrelated(a1,a2,r2):
    
    # CALCULATE: monthly normals and standard errors
    
    x1r = []
    SE1r = []

    for i in range(12):
    
        
        x1a = np.nanmean( a1[a1.index.month==(i+1)] )
        x2a = np.nanmean( a2[a2.index.month==(i+1)] )
        x2r = np.nanmean( r2[r2.index.month==(i+1)] )    
        SE1a = np.nanstd( a1[a1.index.month==(i+1)] ) / np.sqrt( np.isfinite(a1[a1.index.month==(i+1)]).sum() )
        SE2a = np.nanstd( a2[a2.index.month==(i+1)] ) / np.sqrt( np.isfinite(a2[a2.index.month==(i+1)]).sum() )
        SE2r = np.nanstd( r2[r2.index.month==(i+1)] ) / np.sqrt( np.isfinite(r2[r2.index.month==(i+1)]).sum() )
   
        # CASE A: uncorrelated errors
        
        x1r_month = x1a + (x2r - x2a)
        x1r.append(x1r_month)
        SE1r_month = np.sqrt( SE1a**2. + SE2r**2. + SE2a**2. )
        SE1r.append(SE1r_month)
                        
    return x1r, SE1r

def calculate_normals_and_SEs(a1,a2,r2):
    
    # CALCULATE: monthly normals and standard errors
    
    x1r = []
    SE1r = []
    SE12 = []
    n12 = [] 

    for i in range(12):
    
        x1a = np.nanmean( a1[a1.index.month==(i+1)] )
        x2a = np.nanmean( a2[a2.index.month==(i+1)] )
        x2r = np.nanmean( r2[r2.index.month==(i+1)] )    
        SE1a = np.nanstd( a1[a1.index.month==(i+1)] ) / np.sqrt( np.isfinite(a1[a1.index.month==(i+1)]).sum() )
        SE2a = np.nanstd( a2[a2.index.month==(i+1)] ) / np.sqrt( np.isfinite(a2[a2.index.month==(i+1)]).sum() )
        SE2r = np.nanstd( r2[r2.index.month==(i+1)] ) / np.sqrt( np.isfinite(r2[r2.index.month==(i+1)]).sum() )
   
        a12 = a1 - a2
        n12a = np.isfinite(a12[a12.index.month==(i+1)]).sum() 
        x12a = np.nanmean( a12[a12.index.month==(i+1)] )
        x1r_month = x2r + x12a 
        x1r.append(x1r_month)
        SE12a = np.nanstd( a12[a12.index.month==(i+1)] ) / np.sqrt( n12a )
        SE1r_month = np.sqrt( SE2r**2. + SE12a**2. )    
        SE1r.append( SE1r_month )
        SE12.append( SE12a )
        n12.append( n12a )
                        
    return x1r, SE1r, SE12, n12

#==============================================================================
# LOAD: GloSAT absolute temperaturs, lasso, filter (QC) and extract reference station dataframe --> df
#==============================================================================

if load_glosat == True:
           
    print('loading temperatures ...')
                
    df_temp = pd.read_pickle('DATA/df_temp.pkl', compression='bz2') # dataframe of GloSAT absolute temperatures in degrees C
    
    # EXTRACT: test_station (lat,lon)
    
    test_station_lat = df_temp[df_temp['stationcode']==test_station]['stationlat'].iloc[0]
    test_station_lon = df_temp[df_temp['stationcode']==test_station]['stationlon'].iloc[0]    
    
    # LASSO: ref_station list (lat,lon,distance) within a distance test_radius of test_station
        
    station_codes = df_temp['stationcode'].unique()
    station_lats = df_temp.groupby('stationcode')['stationlat'].mean()
    station_lons = df_temp.groupby('stationcode')['stationlon'].mean()
    da = pd.DataFrame({'lat':station_lats, 'lon':station_lons})    
    lasso = find_nearest_lasso( test_station_lat, test_station_lon, da, test_radius)
    lasso['distance'] = [ haversine(test_station_lat, test_station_lon, lasso['lat'].iloc[i], lasso['lon'].iloc[i]) for i in range(len(lasso)) ]

    # CONSTRUCT: dataframe of ref_stations --> df
    
    stationcode_list_lasso = lasso.index.to_list()    
    stationname_list_lasso = [ df_temp[df_temp['stationcode']==stationcode_list_lasso[i]]['stationname'].iloc[0] for i in range(len(stationcode_list_lasso)) ]        

    dates = pd.date_range(start='1700-01-01', end='2021-12-01', freq='MS')
    df = pd.DataFrame(index=dates)
    for i in range(len(stationcode_list_lasso)):                
        dt = df_temp[df_temp['stationcode']==stationcode_list_lasso[i]]
        ts = np.array(dt.groupby('year').mean().iloc[:,0:12]).ravel()                
        if use_fahrenheit == True:             
            ts = centigrade_to_fahrenheit(ts)                    
        t = pd.date_range(start=str(dt.year.iloc[0]), periods=len(ts), freq='MS')
        df[stationcode_list_lasso[i]] = pd.DataFrame({stationname_list_lasso[i]:ts}, index=t) 
        
    # EXTRACT: segment and normal
        
    df_segment = df[ (df.index>=segment_start) & (df.index<=segment_end) ]
    df_normal = df[ (df.index>=normal_start) & (df.index<=normal_end) ]
                    
    # KEEP:   ( stations with > 15 years for all monthly normals in 1961-1990 ) 
    #       & ( ref vs test station overlaps > 5 years for all months in segment ) 
    
    stationcode_list = []
    stationname_list = []
    for station in range(len(stationcode_list_lasso)):
        
        ref_station = stationcode_list_lasso[station]        
        a1 = df_segment[test_station]
        a2 = df_segment[ref_station]
        r2 = df_normal[ref_station]
    
        # COUNT: number of years available to calculate monthly normal and number of years segment overlap
    
        a12_n = []
        r2_n = []
        for i in range(12):            
            n_a = np.isfinite( a1[a1.index.month==(i+1)] + a2[a2.index.month==(i+1)] ).sum()
            n_r = np.isfinite( r2[r2.index.month==(i+1)] ).sum()
            a12_n.append(n_a)
            r2_n.append(n_r)
    
        if ( (np.array(r2_n) > n_baseline_years).sum() == 12 ) & ( (np.array(a12_n) > n_overlap_years).sum() == 12 ):
            stationcode_list.append( stationcode_list_lasso[station] )    
            stationname_list.append( stationname_list_lasso[station] )    

    # DROP: filtered stations from dataframe

    stationcode_list_excluded = list(set(stationcode_list_lasso)-set(stationcode_list))
    stationname_list_excluded = [ df_temp[df_temp['stationcode']==stationcode_list_excluded[i]]['stationname'].iloc[0] for i in range(len(stationcode_list_excluded)) ]                
    lasso_filtered = lasso.T.drop( columns=stationcode_list_excluded ).T
    df_filtered = df.drop( columns=stationcode_list_excluded )

    print('Excluded stations = ', stationcode_list_excluded, stationname_list_excluded)

    # PLOT: distance histogram

    titlestr = 'Reference station distance histogram: bin width = 20 km'
    figstr = 'ref-station-distance-histogram.png'
    
    fig,axs = plt.subplots(figsize=(15,10))
    sns.histplot(lasso_filtered['distance'], stat="probability", binwidth=20, log_scale=False, element="step", fill=True, kde=True)    
    axs.set_xlabel('Distance, km', fontsize=fontsize)
    axs.set_ylabel('Probability', fontsize=fontsize)
    axs.set_title(titlestr, fontsize=fontsize)
    axs.tick_params(labelsize=fontsize)    
    fig.tight_layout()
    plt.savefig(figstr, dpi=300)
    plt.close('all')    
       
#------------------------------------------------------------------------------
# RECALCULATE: test station and neighbouring station mean timeseries in the segment and normal
#------------------------------------------------------------------------------

df_segment = df_filtered[ (df_filtered.index>=segment_start) & (df_filtered.index<=segment_end) ]
df_normal = df_filtered[ (df_filtered.index>=normal_start) & (df_filtered.index<=normal_end) ]

df_test_station = df_filtered[test_station]
df_test_station_segment = df_segment[test_station]
df_test_station_normal = df_normal[test_station]

df_neighbours = df.drop( [test_station], axis=1 )
df_neighbours_segment = df_segment.drop( [test_station], axis=1 )
df_neighbours_normal = df_normal.drop( [test_station], axis=1 )
df_neighbours_mean = df_neighbours.mean(axis=1)
df_neighbours_mean_segment = df_neighbours_segment.mean(axis=1)
df_neighbours_mean_normal = df_neighbours_normal.mean(axis=1)
    
#------------------------------------------------------------------------------
# LOOP: over ref_stations and estimate normal for each station within test_radius
#------------------------------------------------------------------------------

df_errors = pd.DataFrame(columns=[
        'test_station', 
        'ref_station', 
        'error_x1r_CASE_1A',
        'error_x1r_CASE_1B',
        'error_SE1r_CASE_1A',
        'error_SE1r_CASE_1B',
        'median_n12',
        ])
df_errors['test_station'] = len(stationcode_list)*[test_station]    

for station in range(len(stationcode_list)):
    
    ref_station = stationcode_list[station]
    df_ref_station = df[ref_station]
    df_ref_station_segment = df_segment[ref_station]
    df_ref_station_normal = df_normal[ref_station]
    
    # CALCULATE: expected 'true' monthly means from test station
    
    r1 = df_test_station_normal
    x1r_truth = []
    SE1r_truth = []
    for i in range(12):
        
        x1r = np.nanmean( r1[r1.index.month==(i+1)] )    
        x1r_truth.append(x1r)
        SE1r = np.nanstd( r1[r1.index.month==(i+1)] ) / np.sqrt( np.isfinite(r1[r1.index.month==(i+1)]).sum() )   
        SE1r_truth.append(SE1r)
        
    #------------------------------------------------------------------------------
    # MODEL 1: CASE 1A & CASE 1B: single reference: x1=test_station, x2=BHO (referemce_station)
    #------------------------------------------------------------------------------
    
    a1 = df_test_station_segment
    a2 = df_ref_station_segment
    r2 = df_ref_station_normal
    
    x1r_CASE_1A, SE1r_CASE_1A, = calculate_normals_and_SEs_uncorrelated(a1,a2,r2)
    x1r_CASE_1B, SE1r_CASE_1B, SE12, n12 = calculate_normals_and_SEs(a1,a2,r2)
    x1r_CASE_1A_normal = pd.Series(np.tile(x1r_CASE_1A, reps=30), index=r2.index)
    x1r_CASE_1B_normal = pd.Series(np.tile(x1r_CASE_1B, reps=30), index=r2.index)
                      
    print(station, n12)
    
    #------------------------------------------------------------------------------
    # STATISTICS: model versus truth monthly mean errors
    #------------------------------------------------------------------------------
    
    error_x1r_CASE_1A = np.nanmean( np.array(x1r_CASE_1A) - np.array(x1r_truth) )
    error_x1r_CASE_1B = np.nanmean( np.array(x1r_CASE_1B) - np.array(x1r_truth) )
    error_SE1r_CASE_1A = np.nanmean( np.array(SE1r_CASE_1A) - np.array(SE1r_truth) )
    error_SE1r_CASE_1B = np.nanmean( np.array(SE1r_CASE_1B) - np.array(SE1r_truth) )
    
    X_1 = df_test_station.dropna().rolling(nsmooth,center=True).mean()                                               # BHO
    X_1a = X_1[ (X_1.index>=segment_start) & (X_1.index<=segment_end) ]                                     # BHO (segment)
    X_1r_truth = X_1[ (X_1.index>=normal_start) & (X_1.index<=normal_end) ]                                 # BHO (normal) truth
    X_1r_estimate_CASE_1A = X_1r_truth + ( np.nanmean( x1r_CASE_1A_normal ) - np.nanmean( X_1r_truth ) )    # BHO (normal) estimate (CASE_1A)
    X_2_CASE_1A = df_ref_station.dropna().rolling(nsmooth,center=True).mean()                                        # single neighbour
    X_2a_CASE_1A = X_2_CASE_1A[ (X_2_CASE_1A.index>=segment_start) & (X_2_CASE_1A.index<=segment_end) ]     # single neighbour (segment)
    X_2r_CASE_1A = X_2_CASE_1A[ (X_2_CASE_1A.index>=normal_start) & (X_2_CASE_1A.index<=normal_end) ]       # single neighbour (normal)            
    error_CASE_1A = X_1r_estimate_CASE_1A[0] -  X_1r_truth[0]                                               # error relative to expected true normal

    X_1 = df_test_station.dropna().rolling(nsmooth,center=True).mean()                                               # BHO
    X_1a = X_1[ (X_1.index>=segment_start) & (X_1.index<=segment_end) ]                                     # BHO (segment)
    X_1r_truth = X_1[ (X_1.index>=normal_start) & (X_1.index<=normal_end) ]                                 # BHO (normal) truth
    X_1r_estimate_CASE_1B = X_1r_truth + ( np.nanmean( x1r_CASE_1B_normal ) - np.nanmean( X_1r_truth ) )    # BHO (normal) estimate (CASE_1B)
    X_2_CASE_1B = df_ref_station.dropna().rolling(nsmooth,center=True).mean()                                        # single neighbour
    X_2a_CASE_1B = X_2_CASE_1B[ (X_2_CASE_1B.index>=segment_start) & (X_2_CASE_1B.index<=segment_end) ]     # single neighbour (segment)
    X_2r_CASE_1B = X_2_CASE_1B[ (X_2_CASE_1B.index>=normal_start) & (X_2_CASE_1B.index<=normal_end) ]       # single neighbour (normal)            
    error_CASE_1B = X_1r_estimate_CASE_1B[0] -  X_1r_truth[0]                                               # error relative to expected true normal
        
    median_n12 = int(np.round( np.median( n12 ) ))
    
    # STORE: model errors

    df_errors['ref_station'].iloc[station] = ref_station
    df_errors['test_station'].iloc[station] = test_station
    df_errors['error_x1r_CASE_1A'].iloc[station] = error_x1r_CASE_1A
    df_errors['error_x1r_CASE_1B'].iloc[station] = error_x1r_CASE_1B
    df_errors['error_SE1r_CASE_1A'].iloc[station] = error_SE1r_CASE_1A
    df_errors['error_SE1r_CASE_1B'].iloc[station] = error_SE1r_CASE_1B
    df_errors['median_n12'].iloc[station] = median_n12
          
    #==============================================================================
    
    if plot_monthly == True: # one plot per ref_station
        
        # PLOT: monthly normals and standard errors
        
        print('plotting X_{1,r} and SE_{1,r} (CASES 1A and 1B) ...')
            
#        if stationcode_list[station] == 725030:
            
        figstr = 'MODEL-1-monthly-x1r-SE1r' + '-' + ref_station + '(' + stationname_list[station].split('/')[0].replace(' ', '_').lower() + ')'  + '-' + test_station + '(bho)' + '.png'        
        titlestr = r'Model 1B Monthly $\bar{X}_{1,r}$ and $SE_{1,r}$: ' + ref_station + ' (' + stationname_list[station].split('/')[0].replace('_', ' ').title() + ') vs test station = ' + test_station + ' (Blue Hill)'
        
        fig, axs = plt.subplots(2,1, figsize=(15,10))
        axs[0].plot(x1r_truth, marker='s', markersize=20, fillstyle='none', linestyle='none', label=r'Truth: $\bar{X}_{1,r}$')
        axs[0].plot(x1r_CASE_1A, marker='^', markersize=10, alpha=0.5, linestyle='none', label=r'Model 1A: $\bar{X}_{1,r}$ $\mu_{\epsilon}$=' + str(np.round(error_x1r_CASE_1A,2)) + '$^{\circ}$' + temperature_unit) 
        axs[0].plot(x1r_CASE_1B, marker='v', markersize=10, alpha=0.5, linestyle='none', label=r'Model 1B: $\bar{X}_{1,r}$ $\mu_{\epsilon}$=' + str(np.round(error_x1r_CASE_1B,2)) + '$^{\circ}$' + temperature_unit) 
        axs[0].plot(n12, marker='*', markersize=10, alpha=0.5, linestyle='none', label=r'N overlap')        
        axs[0].set_xticks(np.arange(0,12))
        axs[0].set_xticklabels(np.arange(1,13))
        axs[0].tick_params(labelsize=fontsize)    
        axs[0].set_xlabel('Month', fontsize=fontsize)
        axs[0].set_ylabel(r'Mean $\bar{X}_{1,r}$, $^{\circ}$' + temperature_unit, fontsize=fontsize)
        axs[0].legend(loc='upper left', ncol=1, markerscale=1, facecolor='lightgrey', framealpha=0.5, fontsize=12)            
        axs[1].plot(SE1r_truth, marker='s', markersize=20, fillstyle='none', linestyle='none', label=r'Truth: $SE_{1,r}$')
        axs[1].plot(SE1r_CASE_1A, marker='^', markersize=10, alpha=0.5, label=r'Model 1A: $SE_{1,r}$ $\mu_{\epsilon}$=' + str(np.round(error_SE1r_CASE_1A,2)) + '$^{\circ}$' + temperature_unit) 
        axs[1].plot(SE1r_CASE_1B, marker='v', markersize=10, alpha=0.5, label=r'Model 1B: $SE_{1,r}$ $\mu_{\epsilon}$=' + str(np.round(error_SE1r_CASE_1B,2)) + '$^{\circ}$' + temperature_unit)
        axs[1].sharex(axs[0]) 
        axs[1].set_xticks(np.arange(0,12))
        axs[1].set_xticklabels(np.arange(1,13))
        axs[1].tick_params(labelsize=fontsize)    
        axs[1].set_xlabel('Month', fontsize=fontsize)
        axs[1].set_ylabel(r'Standard error $SE_{1,r}$ , $^{\circ}$' + temperature_unit, fontsize=fontsize)
        axs[1].legend(loc='upper left', ncol=1, markerscale=1, facecolor='lightgrey', framealpha=0.5, fontsize=12)            
        if use_fahrenheit == True:
            axs[0].set_ylim(0,80)
            axs[1].set_ylim(0,2)
        else:
            ax.set_xlim(-20,40)
            ax.set_ylim(0,0.05)    
        plt.suptitle(titlestr)        
        fig.tight_layout()
        plt.savefig(figstr, dpi=300)
        plt.close('all')
    
    if plot_fit == True: # one plot per ref_station
        
        print('plotting fit X1r ...')
        print('plotting MODEL 1 (nearest neighbour) X_{1,r} ...')            
           
        figstr = 'MODEL-1-fit' + '-' + ref_station + '(' + stationname_list[station].split('/')[0].replace(' ', '_').lower() + ')' + '-' + test_station + '(bho)' + '.png'        
        titlestr = r'Model 1B fit: ' + ref_station + ' (' + stationname_list[station].split('/')[0].replace('_', ' ').title() + ') vs test station = ' + test_station + ' (Blue Hill)'
           
        fig, axs = plt.subplots(figsize=(15,10))
        axs.plot(X_1.index, X_1, marker='.', color='lightgrey', linestyle='None', alpha=1, label='$T_{g}$ ' + test_station)
        axs.plot(X_2a_CASE_1A.index, X_2a_CASE_1A, marker='.', linestyle='None',color='lightblue', alpha=1, label='$X_{2,a}$ ' + ref_station + ' (segment)')
        axs.plot(X_2r_CASE_1A.index, X_2r_CASE_1A, marker='.', linestyle='None',color='pink', alpha=1, label='$X_{2,r}$ ' + ref_station + ' (1961-1990)')
        axs.plot(X_1a.index, X_1a, marker='.', linestyle='None',color='blue', alpha=1, label='$X_{1,a}$ ' + test_station + ' (segment)')
        axs.plot(X_1r_truth.index, X_1r_truth, marker='.', linestyle='None',color='red', alpha=1, label='$X_{1,r}$ ' + test_station + ' (1961-1990) truth')
        axs.plot(X_1r_truth.index, X_1r_estimate_CASE_1A, marker='.', linestyle='None',color='k', alpha=1, label='$X_{1,r}$ ' + test_station + ' (1961-1990) estimate: $\epsilon$=' + str(np.round(error_CASE_1A,2)) + '$^{\circ}$' + temperature_unit)   
        axs.plot(X_2a_CASE_1A.index, len(X_2a_CASE_1A)*[ np.nanmean( X_2a_CASE_1A ) ], ls='--', lw=2, color='lightblue', alpha=1)            
        axs.plot(X_2r_CASE_1A.index, len(X_2r_CASE_1A)*[ np.nanmean( X_2r_CASE_1A ) ], ls='--', lw=2, color='pink', alpha=1)            
        axs.plot(X_1a.index, len(X_1a)*[ np.nanmean( X_1a ) ], ls='--', lw=2, color='blue', alpha=1)            
        axs.plot(X_1r_truth.index, len(X_1r_truth)*[ np.nanmean( X_1r_truth ) ], ls='--', lw=2, color='red', alpha=1)            
        axs.plot(X_1r_truth.index, len(X_1r_estimate_CASE_1A)*[ np.nanmean( x1r_CASE_1A_normal ) ], ls='--', lw=2, color='k', alpha=1)           
        axs.legend(loc='upper left', ncol=1, markerscale=1, facecolor='lightgrey', framealpha=0.5, fontsize=12)        
        axs.set_xlabel('Year', fontsize=fontsize)
        axs.set_ylabel(r'Average absolute temperature (5yr MA), $^{\circ}$' + temperature_unit, fontsize=fontsize)
        axs.set_title(titlestr, fontsize=fontsize)
        axs.tick_params(labelsize=fontsize)   
        if use_fahrenheit == True:
            axs.set_ylim(35,60)
        else:
            ax.set_ylim(-20,40)
        fig.tight_layout()
        plt.savefig(figstr, dpi=300)
        plt.close('all')

# SAVE: model errors array
        
df_errors.to_csv('MODEL-1-errors.csv')
    
#------------------------------------------------------------------------------
print('** END')




