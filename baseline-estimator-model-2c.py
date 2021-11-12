#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#------------------------------------------------------------------------------
# PROGRAM: baseline-estimator-model-2c.py
#------------------------------------------------------------------------------
# Verion 0.3
# 12 July, 2021
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
plot_methodology = True

use_fahrenheit = True
if use_fahrenheit: 
    temperature_unit = 'F'
else:
    temperature_unit = 'C'

#segment_start = pd.to_datetime('1851-01-01')
segment_end = pd.to_datetime('1920-12-01')
normal_start = pd.to_datetime('1961-01-01')
normal_end = pd.to_datetime('1990-12-01')

test_station = '744920'    # BHO
#test_station = '103810'    # Berlin-Dahlem
#test_station = '109620'    # Hohenpeissenberg
#test_station = '024581'    # Uppsala
#test_station = '113200'    # Innsbruck
#test_station = '037401'    # CET

test_radius = 312           # km
n_baseline_years = 15       # Minimum number of years in baseline for calculation of normals 
#n_overlap_years = 25        # Minimum number of years overlap in the segment epoch

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
                        
    return x1r, SE1r, SE12, n12, SE1a, SE2a, SE2r

#==============================================================================
# LOAD: GloSAT absolute temperaturs, lasso, filter (QC) and extract reference station dataframe --> df
#==============================================================================

if load_glosat == True:
           
    print('loading temperatures ...')
                
    df_temp = pd.read_pickle('DATA/df_temp.pkl', compression='bz2') # dataframe of GloSAT absolute temperatures in degrees C
        
    # FIX: set pandas >1678 to fix datetime calendar limit

    df_temp = df_temp[ df_temp.year>1678 ]
    
    # EXTRACT: test_station (lat,lon)
    
    test_station_lat = df_temp[df_temp['stationcode']==test_station]['stationlat'].iloc[0]
    test_station_lon = df_temp[df_temp['stationcode']==test_station]['stationlon'].iloc[0]    
    test_station_name = df_temp[df_temp['stationcode']==test_station]['stationname'].iloc[0]    
    
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
        
    df_normal = df[ (df.index>=normal_start) & (df.index<=normal_end) ]
        
    # KEEP: ( stations with > 15 years for all monthly normals in 1961-1990 ) 
    
    stationcode_list = []
    stationname_list = []
    for station in range(len(stationcode_list_lasso)):
        
        ref_station = stationcode_list_lasso[station]        
        r2 = df_normal[ref_station]    
        r2_n = []
        for i in range(12):            
            n_r = np.isfinite( r2[r2.index.month==(i+1)] ).sum()
            r2_n.append(n_r)
    
        if (np.array(r2_n) > n_baseline_years).sum() == 12:
            stationcode_list.append( stationcode_list_lasso[station] )    
            stationname_list.append( stationname_list_lasso[station] )    

    # DROP: filtered stations from dataframe

    stationcode_list_excluded = list(set(stationcode_list_lasso)-set(stationcode_list))
    stationname_list_excluded = [ df_temp[df_temp['stationcode']==stationcode_list_excluded[i]]['stationname'].iloc[0] for i in range(len(stationcode_list_excluded)) ]                
    lasso_filtered = lasso.T.drop( columns=stationcode_list_excluded ).T
    df_filtered = df.drop( columns=stationcode_list_excluded )
    df_filtered_baseline = df_filtered[ (df_filtered.index >= normal_start) & (df_filtered.index <= normal_end) ]
    df_filtered_external = df_filtered[ df_filtered.index < segment_end ]

    mask = np.isfinite( df_filtered_external )      
    n_stations = mask.sum(axis=1)   
    n_months = mask.sum(axis=0)       
    mask_n_months = n_months < (100*12)
    mask_n_stations = n_stations < 10    
    df_filtered_external_n_months = df_filtered_external.drop( columns=mask_n_months.index[mask_n_months] )
    df_filtered_external_n_months_n_stations = df_filtered_external_n_months.T.drop( columns=mask_n_stations.index[mask_n_stations] ).T.dropna()

    stationcode_list = df_filtered_external_n_months_n_stations.columns.values
    stationcode_list_excluded = list(set(stationcode_list_lasso)-set(stationcode_list))
    stationname_list_excluded = [ df_temp[df_temp['stationcode']==stationcode_list_excluded[i]]['stationname'].iloc[0] for i in range(len(stationcode_list_excluded)) ]                
    lasso_filtered = lasso.T.drop( columns=stationcode_list_excluded ).T
    df_filtered = df.drop( columns=stationcode_list_excluded )
    
    print('Excluded stations = ', stationcode_list_excluded, stationname_list_excluded)

    df_excluded_stations = pd.DataFrame({'stationcode_list_excluded':stationcode_list_excluded, 'stationname_list_excluded':stationname_list_excluded})
    df_excluded_stations.to_csv('MODEL-2C-excluded-stations' + '-' + test_station + '(' + test_station_name.replace(' ','_') + ')' + '-' + str(test_radius) + 'km' + '.csv')
                
#------------------------------------------------------------------------------
# RECALCULATE: test station and neighbouring station mean timeseries in the segment and normal
#------------------------------------------------------------------------------

df_normal = df_filtered[ (df_filtered.index >= normal_start) & (df_filtered.index <= normal_end) ]
#df_segment = df_filtered[ (df_filtered.index >= segment_start) & (df_filtered.index <= segment_end) ]
df_segment = df_filtered[ (df_filtered.index < segment_end) ].dropna()
segment_start = df_segment.index[0]
#segment_end = df_segment.index[-1]
df_test_station = df_filtered[test_station]
df_test_station_segment = df_test_station[ (df_test_station.index >= segment_start) & (df_test_station.index <= segment_end) ]
df_test_station_normal = df_test_station[ (df_test_station.index >= normal_start) & (df_test_station.index <= normal_end) ]

df_neighbours = df_filtered.drop( [test_station], axis=1 ).dropna()
df_neighbours_segment = df_neighbours[ (df_neighbours.index >= segment_start) & (df_neighbours.index <= segment_end) ]
df_neighbours_normal = df_neighbours[ (df_neighbours.index >= normal_start) & (df_neighbours.index <= normal_end) ]
df_neighbours_mean = df_neighbours.mean(axis=1)
df_neighbours_mean_segment = df_neighbours_segment.mean(axis=1)
df_neighbours_mean_normal = df_neighbours_normal.mean(axis=1)

#------------------------------------------------------------------------------
# CALCULATE: Model 2C predictor using mean of core stations and test_station segment
#------------------------------------------------------------------------------

df_errors = pd.DataFrame(columns=[
        'error_x1r',
        'error_SE1r',
        'error_SE12',
        ], 
    index=[test_station])
        
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
# MODEL 2C: x1=mean of neighbour sub-ensemble, x2=test_station
#------------------------------------------------------------------------------
    
a1 = df_test_station_segment
a2 = df_neighbours_mean_segment
r2 = df_neighbours_mean_normal    

x1r, SE1r, SE12, n12, SE1a, SE2a, SE2r = calculate_normals_and_SEs(a1,a2,r2)
x1r_normal = pd.Series(np.tile(x1r, reps=30), index=pd.date_range(start='1961-01-01', periods=360, freq='MS'))
    
#------------------------------------------------------------------------------
# STATISTICS: model versus truth monthly mean errors
#------------------------------------------------------------------------------
               
error_x1r = np.nanmean( np.array(x1r) - np.array(x1r_truth) )
error_SE1r = np.nanmean( np.array(SE1r) - np.array(SE1r_truth) )
error_SE12 = np.nanmean( np.array(SE12) )

X_1 = df_test_station.rolling(nsmooth,center=True).mean()                                           # test_station
#X_1a = X_1[ (X_1.index>=segment_start) & (X_1.index<=segment_end) ]                                # test_station (segment)
X_1a = X_1[ (X_1.index < segment_end) ]                                                             # test_station (segment)
X_1r_truth = X_1[ (X_1.index >= normal_start) & (X_1.index <= normal_end) ]                         # test_station (normal) truth
X_1r_estimate = X_1r_truth + ( np.nanmean( x1r_normal ) - np.nanmean( X_1r_truth ) )                # test_station (normal) estimate (Model 2C)
error = X_1r_estimate[0] -  X_1r_truth[0]                                                    	    # error relative to expected true normal

X_2 = df_neighbours_mean.rolling(nsmooth,center=True).mean()                                        # Neighbours mean
#X_2a = X_2[ (X_2.index >= segment_start) & (X_2.index <= segment_end) ]                            # Neighbours mean (segment)
X_2a = X_2[ (X_2.index < segment_end) ].dropna()                                                    # Neighbours mean (segment)
X_2r = X_2[ (X_2.index >= normal_start) & (X_2.index <= normal_end) ]                               # Neighbours mean (normal)

X = df_neighbours.rolling(nsmooth,center=True).mean()                                               # All neighbours
#X_segment = X[ (X.index>=segment_start) & (X.index<=segment_end) ]                                 # All neighbours (segment)
X_segment = X[ (X.index < segment_end) ].dropna()                                                   # All neighbours (segment)
X_normal = X[ (X.index >= normal_start) & (X.index <= normal_end) ]                                 # All neighbours (normal)

# SAVE: model errors array
    
df_errors['error_x1r'].iloc[0] = error_x1r
df_errors['error_SE1r'].iloc[0] = error_SE1r
df_errors['error_SE12'].iloc[0] = error_SE12
df_errors.to_csv('MODEL-2C-errors' + '-' + test_station + '(' + test_station_name.replace(' ','_') + ')' + '-' + str(test_radius) + 'km' + '.csv')

#==============================================================================
            
# PLOT: monthly normals and standard errors
        
print('plotting X_{1,r} and SE_{1,r} ...')
            
figstr = 'MODEL-2C-monthly-x1r-SE1r' + '-' + test_station + '(' + test_station_name.replace(' ','_') + ')' + '-' + str(test_radius) + 'km' + '.png'
titlestr = r'Model 2C monthly $\bar{X}_{1,r}$ and $SE_{1,r}$: ' + str(df_neighbours.shape[1]) + ' neighbours < ' + str(test_radius) + ' km of the test station = ' + test_station + ' (' + test_station_name.replace('_',' ').title() + ')'
        
fig, axs = plt.subplots(2,1, figsize=(15,10))
axs[0].plot(x1r_truth, marker='s', markersize=20, fillstyle='none', linestyle='none', label=r'Truth: $\bar{X}_{1,r}$')
axs[0].plot(x1r, marker='o', markersize=10, alpha=0.5, linestyle='none', label=r'Model 2C: $\bar{X}_{1,r}$ $\mu_{\epsilon}$=' + str(np.round(error_x1r,2)) + '$^{\circ}$' + temperature_unit) 
axs[0].set_xticks(np.arange(0,12))
axs[0].set_xticklabels(np.arange(1,13))
axs[0].tick_params(labelsize=fontsize)    
axs[0].set_xlabel('Month', fontsize=fontsize)
axs[0].set_ylabel(r'Mean $\bar{X}_{1,r}$, $^{\circ}$' + temperature_unit, fontsize=fontsize)
axs[0].legend(loc='upper left', ncol=1, markerscale=1, facecolor='lightgrey', framealpha=0.5, fontsize=12)            
axs[1].plot(SE1r_truth, marker='s', markersize=20, fillstyle='none', linestyle='none', label=r'Truth: $SE_{1,r}$')
axs[1].plot(SE1r, marker='o', markersize=10, alpha=0.5, label=r'Model 2C: $SE_{1,r}$ $\mu_{\epsilon}$=' + str(np.round(error_SE1r,2)) + '$^{\circ}$' + temperature_unit)
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
plt.suptitle(titlestr, fontsize=fontsize)        
fig.tight_layout()
plt.savefig(figstr, dpi=300)
plt.close('all')
    
print('plotting neighbours X_{1,r} ...')
                
figstr = 'MODEL-2C-fit' + '-' + test_station + '(' + test_station_name.replace(' ','_') + ')' + '-' + str(test_radius) + 'km' + '.png'
titlestr = r'Model 2C fit: $X_{1,r}$: ' + str(df_neighbours.shape[1]) + ' neighbours < ' + str(test_radius) + ' km of the test station = ' + test_station.title() + ' (' + test_station_name.replace('_',' ').title() + ')'
        
fig, axs = plt.subplots(figsize=(15,10))
axs.plot(X_segment.index, X_segment, marker='.', ls='none', color='lightblue', alpha=1)
axs.plot(X_normal.index, X_normal, marker='.', ls='none', color='pink', alpha=1)
axs.plot(X_1.index, X_1, marker='.', ls='none', color='lightgrey', alpha=1, label='$T_{g}$ ' + test_station)
axs.plot(X_2r.index, X_2r, marker='.', ls='none', color='purple', alpha=1, label='$X_{2,r}$ neighbours mean (1961-1990)')
axs.plot(X_2a.index[X_2a.index>X_segment.index[0]], X_2a[X_2a.index>X_segment.index[0]], marker='.', ls='none', color='teal', alpha=1, label='$X_{2,a}$ neighbours mean (segment)')
axs.plot(X_1a.index[X_1a.index>X_segment.index[0]], X_1a[X_1a.index>X_segment.index[0]], marker='.', ls='none', color='blue', alpha=1, label='$X_{1,a}$ ' + test_station + ' (segment)')
axs.plot(X_1r_truth.index, X_1r_truth, marker='.', ls='none', color='red', alpha=1, label='$X_{1,r}$ ' + test_station + ' (1961-1990) truth')    
axs.plot(X_1r_estimate.index, X_1r_estimate, marker='.', ls='none', color='k', alpha=1, label='$X_{1,r}$ ' + test_station + ' (1961-1990) estimate: $\epsilon$=' + str(np.round(error,2)) + '$^{\circ}$' + temperature_unit)    
axs.plot(X_2r.index, len(X_2r)*[ np.nanmean( X_2r ) ], ls='--', lw=2, color='purple', alpha=1)            
axs.plot(X_2a.index[X_2a.index>X_segment.index[0]], len(X_2a[X_2a.index>X_segment.index[0]])*[ np.nanmean( X_2a ) ], ls='--', lw=2, color='teal', alpha=1)            
axs.plot(X_1a.index[X_1a.index>X_segment.index[0]], len(X_1a[X_1a.index>X_segment.index[0]])*[ np.nanmean( X_1a ) ], ls='--', lw=2, color='blue', alpha=1)            
axs.plot(X_1r_truth.index, len(X_1r_truth)*[ np.nanmean( X_1r_truth ) ], ls='--', lw=2, color='red', alpha=1)            
axs.plot(x1r_normal.index, len(x1r_normal)*[ np.nanmean( x1r_normal ) ], ls='--', lw=2, color='k', alpha=1)           
axs.legend(loc='upper left', ncol=1, markerscale=1, facecolor='lightgrey', framealpha=0.5, fontsize=12)        
axs.set_xlabel('Year', fontsize=fontsize)
axs.set_ylabel(r'Average absolute temperature (5yr MA), $^{\circ}$' + temperature_unit, fontsize=fontsize)
axs.set_title(titlestr, fontsize=fontsize)
axs.tick_params(labelsize=fontsize)   
if use_fahrenheit == True:
    axs.set_xlim(pd.to_datetime('1700-01-01'), pd.to_datetime('2021-12-01'))
    axs.set_ylim(35,60)
else:
    axs.set_xlim(pd.to_datetime('1700-01-01'), pd.to_datetime('2021-12-01'))
    axs.set_ylim(-20,40)
fig.tight_layout()
plt.savefig(figstr, dpi=300)
plt.close('all')            
   
#==============================================================================

if plot_methodology == True:            
            
    print('plotting the methodology ...')
                    
    # Segment X_{1,a}
    
    figstr = 'methodology-x1a.png'
    titlestr = r'Methodology: thought experiment - let us take a short (real) segment $X_{1,a}$ from HadCET'
            
    fig, axs = plt.subplots(figsize=(15,10))
    axs.plot(X_1a.index[X_1a.index>X_segment.index[0]], X_1a[X_1a.index>X_segment.index[0]], marker='.', ls='none', color='blue', alpha=1, label='$X_{1,a}$ ' + test_station + ' (segment)')
    axs.legend(loc='upper left', ncol=1, markerscale=1, facecolor='lightgrey', framealpha=0.5, fontsize=12)        
    axs.set_xlabel('Year', fontsize=fontsize)
    axs.set_ylabel(r'Average absolute temperature (5yr MA), $^{\circ}$' + temperature_unit, fontsize=fontsize)
    axs.set_title(titlestr, fontsize=fontsize)
    axs.tick_params(labelsize=fontsize)   
    if use_fahrenheit == True:
        axs.set_xlim(pd.to_datetime('1700-01-01'), pd.to_datetime('2021-12-01'))
        axs.set_ylim(35,60)
    else:
        ax.set_ylim(-20,40)
    fig.tight_layout()
    plt.savefig(figstr, dpi=300)
    plt.close('all')           

    # Segment X_{1,a} and desired normal X_{1,r}

    figstr = 'methodology-x1a-x1r.png'
    titlestr = r'Methodology: thought experiment - let us say we want to estimate its unknown monthly normals in $X_{1,r}$'
            
    fig, axs = plt.subplots(figsize=(15,10))
    axs.plot(X_1a.index[X_1a.index>X_segment.index[0]], X_1a[X_1a.index>X_segment.index[0]], marker='.', ls='none', color='blue', alpha=1, label='$X_{1,a}$ ' + test_station + ' (segment)')
    axs.plot(X_1r_truth.index, X_1r_truth, marker='.', ls='none', color='red', alpha=1, label='$X_{1,r}$ ' + test_station + ' (1961-1990) desired normal')    
    axs.legend(loc='upper left', ncol=1, markerscale=1, facecolor='lightgrey', framealpha=0.5, fontsize=12)        
    axs.set_xlabel('Year', fontsize=fontsize)
    axs.set_ylabel(r'Average absolute temperature (5yr MA), $^{\circ}$' + temperature_unit, fontsize=fontsize)
    axs.set_title(titlestr, fontsize=fontsize)
    axs.tick_params(labelsize=fontsize)   
    if use_fahrenheit == True:
        axs.set_xlim(pd.to_datetime('1700-01-01'), pd.to_datetime('2021-12-01'))
        axs.set_ylim(35,60)
    else:
        ax.set_ylim(-20,40)
    fig.tight_layout()
    plt.savefig(figstr, dpi=300)
    plt.close('all')    

    # Segment X_{1,a}, reference segment X_{2,a}, reference normal X_{2,r}

    figstr = 'methodology-x1a-x2a-x2r.png'
    titlestr = r'Methodology: we have a neighbouring reference series with co-located observations in the segment $X_{2,a}$ and baseline $X_{2,r}$'
            
    fig, axs = plt.subplots(figsize=(15,10))
    axs.plot(X.index, X, marker='.', ls='none', color='lightgrey', alpha=1)
    axs.plot(X_1a.index[X_1a.index>X_segment.index[0]], X_1a[X_1a.index>X_segment.index[0]], marker='.', ls='none', color='blue', alpha=1, label='$X_{1,a}$ ' + test_station + ' (segment)')
    axs.plot(X_2a.index[X_2a.index>X_segment.index[0]], X_2a[X_2a.index>X_segment.index[0]], marker='.', ls='none', color='teal', alpha=1, label='$X_{2,a}$ neighbour (segment)')
    axs.plot(X_2r.index, X_2r, marker='.', ls='none', color='purple', alpha=1, label='$X_{2,r}$ neighbour (1961-1990)')
    axs.legend(loc='upper left', ncol=1, markerscale=1, facecolor='lightgrey', framealpha=0.5, fontsize=12)        
    axs.set_xlabel('Year', fontsize=fontsize)
    axs.set_ylabel(r'Average absolute temperature (5yr MA), $^{\circ}$' + temperature_unit, fontsize=fontsize)
    axs.set_title(titlestr, fontsize=fontsize)
    axs.tick_params(labelsize=fontsize)   
    if use_fahrenheit == True:
        axs.set_xlim(pd.to_datetime('1700-01-01'), pd.to_datetime('2021-12-01'))
        axs.set_ylim(35,60)
    else:
        ax.set_ylim(-20,40)
    fig.tight_layout()
    plt.savefig(figstr, dpi=300)
    plt.close('all')           

    # Segment X_{1,a}, Xreference segment X_{2,a}, reference normal X_{2,r} and desired normal X_{1,r}

    figstr = 'methodology-x1a-x2a-x2r-x1r.png'
    titlestr = r'Methodology: can we make a model to estimate $X_{1,r}$ from $X_{1,a}$, $X_{2,a}$ and $X_{2,r}$ ?'
            
    fig, axs = plt.subplots(figsize=(15,10))
    axs.plot(X_1a.index[X_1a.index>X_segment.index[0]], X_1a[X_1a.index>X_segment.index[0]], marker='.', ls='none', color='blue', alpha=1, label='$X_{1,a}$ ' + test_station + ' (segment)')
    axs.plot(X_2a.index[X_2a.index>X_segment.index[0]], X_2a[X_2a.index>X_segment.index[0]], marker='.', ls='none', color='teal', alpha=1, label='$X_{2,a}$ neighbour (segment)')
    axs.plot(X_2r.index, X_2r, marker='.', ls='none', color='purple', alpha=1, label='$X_{2,r}$ neighbour (1961-1990)')
    axs.plot(X_1r_truth.index, X_1r_truth, marker='.', ls='none', color='red', alpha=1, label='$X_{1,r}$ ' + test_station + ' (1961-1990) desired normal')    
    axs.legend(loc='upper left', ncol=1, markerscale=1, facecolor='lightgrey', framealpha=0.5, fontsize=12)        
    axs.set_xlabel('Year', fontsize=fontsize)
    axs.set_ylabel(r'Average absolute temperature (5yr MA), $^{\circ}$' + temperature_unit, fontsize=fontsize)
    axs.set_title(titlestr, fontsize=fontsize)
    axs.tick_params(labelsize=fontsize)   
    if use_fahrenheit == True:
        axs.set_xlim(pd.to_datetime('1700-01-01'), pd.to_datetime('2021-12-01'))
        axs.set_ylim(35,60)
    else:
        ax.set_ylim(-20,40)
    fig.tight_layout()
    plt.savefig(figstr, dpi=300)
    plt.close('all')         
    
    # Segment X_{1,a}, reference segment X_{2,a}, reference normal X_{2,r}

    figstr = 'methodology-x1a-x2a-x2r-x1r-abstracted-means.png'
    titlestr = r'Methodology: we abstract and develop a model based on the mean values of the timeseries levels'
            
    fig, axs = plt.subplots(figsize=(15,10))
    axs.plot(X_1a.index[X_1a.index>X_segment.index[0]], len(X_1a[X_1a.index>X_segment.index[0]])*[ np.nanmean( X_1a ) ], ls='--', lw=2, color='blue', alpha=1, label='$X_{1,a}$ ' + test_station + ' (segment)')            
    axs.plot(X_2a.index[X_2a.index>X_segment.index[0]], len(X_2a[X_2a.index>X_segment.index[0]])*[ np.nanmean( X_2a ) ], ls='--', lw=2, color='teal', alpha=1, label='$X_{2,a}$ neighbour (segment)')            
    axs.plot(X_2r.index, len(X_2r)*[ np.nanmean( X_2r ) ], ls='--', lw=2, color='purple', alpha=1, label='$X_{2,r}$ neighbour (1961-1990)')            
    axs.plot(x1r_normal.index, len(x1r_normal)*[ np.nanmean( x1r_normal ) ], ls='--', lw=2, color='k', alpha=1, label='$X_{1,r}$ ' + test_station + ' (1961-1990) estimate: $\epsilon$=' + str(np.round(error,2)) + '$^{\circ}$' + temperature_unit)    
    axs.plot(X_1r_truth.index, len(X_1r_truth)*[ np.nanmean( X_1r_truth ) ], ls='--', lw=2, color='red', alpha=1, label='$X_{1,r}$ ' + test_station + ' (1961-1990) desired normal')                
    axs.legend(loc='upper left', ncol=1, markerscale=1, facecolor='lightgrey', framealpha=0.5, fontsize=12)        
    axs.set_xlabel('Year', fontsize=fontsize)
    axs.set_ylabel(r'Average absolute temperature (5yr MA), $^{\circ}$' + temperature_unit, fontsize=fontsize)
    axs.set_title(titlestr, fontsize=fontsize)
    axs.tick_params(labelsize=fontsize)   
    if use_fahrenheit == True:
        axs.set_xlim(pd.to_datetime('1700-01-01'), pd.to_datetime('2021-12-01'))
        axs.set_ylim(35,60)
    else:
        ax.set_ylim(-20,40)
    fig.tight_layout()
    plt.savefig(figstr, dpi=300)
    plt.close('all')         

    # Segment X_{1,a}, reference segment X_{2,a}, reference normal X_{2,r}

    figstr = 'methodology-x1a-x2a-x2r-x1r-abstracted-means-sterr.png'
    titlestr = r'Methodology: we model out correlations and combine standard errors in quadrature to calculate the uncertainty on the model estimated mean'
            
    fig, axs = plt.subplots(figsize=(15,10))
    axs.plot(X_1a.index[X_1a.index>X_segment.index[0]], len(X_1a[X_1a.index>X_segment.index[0]])*[ np.nanmean( X_1a ) ], ls='--', lw=2, color='blue', alpha=1, label='$X_{1,a}$ ' + test_station + ' (segment)')            
    axs.plot(X_2a.index[X_2a.index>X_segment.index[0]], len(X_2a[X_2a.index>X_segment.index[0]])*[ np.nanmean( X_2a ) ], ls='--', lw=2, color='teal', alpha=1, label='$X_{2,a}$ neighbour (segment)')            
    axs.plot(X_2r.index, len(X_2r)*[ np.nanmean( X_2r ) ], ls='--', lw=2, color='purple', alpha=1, label='$X_{2,r}$ neighbour (1961-1990)')            
    axs.plot(x1r_normal.index, len(x1r_normal)*[ np.nanmean( x1r_normal ) ], ls='--', lw=2, color='k', alpha=1, label='$X_{1,r}$ ' + test_station + ' (1961-1990) estimate: $\epsilon$=' + str(np.round(error,2)) + '$^{\circ}$' + temperature_unit)    
 
    axs.fill_between(X_1a.index[X_1a.index>X_segment.index[0]], len(X_1a[X_1a.index>X_segment.index[0]])*[ np.nanmean( X_1a ) - SE1a ], len(X_1a[X_1a.index>X_segment.index[0]])*[ np.nanmean( X_1a ) + SE1a ], color='blue', alpha=0.2, label='$SE_{1,a}$ ' + test_station + ' (segment)')            
    axs.fill_between(X_2a.index[X_2a.index>X_segment.index[0]], len(X_2a[X_2a.index>X_segment.index[0]])*[ np.nanmean( X_2a ) - SE2a ], len(X_2a[X_2a.index>X_segment.index[0]])*[ np.nanmean( X_2a ) + SE2a ], color='teal', alpha=0.2, label='$SE_{2,a}$ neighbour (segment)')            
    axs.fill_between(X_2r.index[X_2r.index>X_segment.index[0]], len(X_2r[X_2r.index>X_segment.index[0]])*[ np.nanmean( X_2r ) - SE2r ], len(X_2r[X_2r.index>X_segment.index[0]])*[ np.nanmean( X_2r ) + SE2r ], color='purple', alpha=0.2, label='$SE_{2,r}$ neighbour (1961-1990)')            
    axs.fill_between(x1r_normal.index[x1r_normal.index>X_segment.index[0]], len(x1r_normal[x1r_normal.index>X_segment.index[0]])*[ np.nanmean( x1r_normal ) - np.nanmean( np.array(SE1r) ) ], len(x1r_normal[x1r_normal.index>X_segment.index[0]])*[ np.nanmean( x1r_normal ) + np.nanmean( np.array(SE1r) ) ], color='black', alpha=0.2, label='$SE_{1,r}$ ' + test_station + ' (1961-1990) estimate')           

    axs.plot(X_1r_truth.index, len(X_1r_truth)*[ np.nanmean( X_1r_truth ) ], ls='--', lw=2, color='red', alpha=1, label='$X_{1,r}$ ' + test_station + ' (1961-1990) desired normal')                
        
    axs.legend(loc='upper left', ncol=1, markerscale=1, facecolor='lightgrey', framealpha=0.5, fontsize=12)        
    axs.set_xlabel('Year', fontsize=fontsize)
    axs.set_ylabel(r'Average absolute temperature (5yr MA), $^{\circ}$' + temperature_unit, fontsize=fontsize)
    axs.set_title(titlestr, fontsize=fontsize)
    axs.tick_params(labelsize=fontsize)   
    if use_fahrenheit == True:
        axs.set_xlim(pd.to_datetime('1700-01-01'), pd.to_datetime('2021-12-01'))
        axs.set_ylim(35,60)
    else:
        ax.set_ylim(-20,40)
    fig.tight_layout()
    plt.savefig(figstr, dpi=300)
    plt.close('all')         
            
#------------------------------------------------------------------------------
print('** END')
