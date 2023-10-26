### FINAL Working global code ####


import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import iris
from iris import plot
import cf_units
import cftime
import cartopy
import cartopy.crs as ccrs
import geopandas as gpd
import warnings
import pandas as pd
import pickle
from copy import deepcopy
from scipy import stats
from scipy.signal import periodogram
from scipy.stats import pearsonr
import random

SEED = int.from_bytes("They're taking the hobbits to Isengard!".encode('utf-8'), byteorder='big')

def constrain_time(cube, date1, date2):
    if not any([coord.long_name == 'year' for coord in cube.coords()]):
        iris.coord_categorisation.add_year(cube, 'time', name='year')
    return cube.extract(iris.Constraint(time=lambda cell: date1 <= cell.point.year <= date2))

def to_timeseries(cube):
    coords = ('longitude', 'latitude')
    for coord in coords:
        if not cube.coord(coord).has_bounds():
            cube.coord(coord).guess_bounds()
    weights = iris.analysis.cartography.area_weights(cube)
    return cube.collapsed(coords, iris.analysis.SUM, weights = weights)/10**12

def to_monthly_mean(cube):
    if not any([coord.long_name == 'month' for coord in cube.coords()]):
        iris.coord_categorisation.add_month(cube, 'time', name='month')
    return cube.aggregated_by('month', iris.analysis.MEAN)

def pixel_mean(cube):
    return cube.collapsed(('time'), iris.analysis.MEAN)

def to_monthly_max(cube):
    if not any([coord.long_name == 'month' for coord in cube.coords()]):
        iris.coord_categorisation.add_month(cube, 'time', name='month')
    return cube.aggregated_by('month', iris.analysis.MAX)

def to_annual(cube):
    if not any([coord.long_name == 'year' for coord in cube.coords()]):
        iris.coord_categorisation.add_year(cube, 'time', name='year')
    return cube.aggregated_by(['year'], iris.analysis.SUM)

def add_year(cube):
    if not any([coord.long_name == 'year' for coord in cube.coords()]):
        iris.coord_categorisation.add_year(cube, 'time', name='year')
    return cube

def to_df(cube, var_name='month', val_name='Burned Area (in Mha)'):
    cube_data = cube.data.compressed().reshape(cube.shape[0], -1)
    cube_df = pd.DataFrame(data=cube_data).transpose().melt()
    cube_df.rename(columns = {'variable':var_name, 'value':val_name}, inplace = True)
    return cube_df

def to_percentiles(cube, percentiles):
    coords = ('longitude', 'latitude')
    for coord in coords:
        if not cube.coord(coord).has_bounds():
            cube.coord(coord).guess_bounds()
    return cube.collapsed(coords, iris.analysis.PERCENTILE, percent=percentiles)

def from_perc_to_mha(cube):
    weights = iris.analysis.cartography.area_weights(cube)
    return iris.analysis.maths.multiply(cube, weights)/10**12

def mask_below_percentile(cube, p=0.99):
    pXX = cube.collapsed(('longitude', 'latitude', 'time'), iris.analysis.PERCENTILE, percent=p)
    mask = xr.where(cube.data < pXX.data, True, False)
    return iris.util.mask_cube(cube, mask)

def count_per_time(masked_cube):
    counts = masked_cube.collapsed(('longitude', 'latitude'), iris.analysis.COUNT, function=lambda values: values > 0)
    return counts

def to_Z_score(cube):
    anomaly = cube-cube.collapsed('time', iris.analysis.MEAN)
    return anomaly/np.std(cube.data)

def df_to_annual(df):
    # 1. Convert the index to year to_period('Y'), which turns it into a pandas.PeriodIndex.
    # 2. Get the year as int (.year, seaborn doesn't know how to handle a pandas.PeriodIndex)
    # 3. Group by year .groupby()
    # 4. Take the sum of each group .sum()
    # 5. NaNs are treated as 0s, and a series of all NaNs would equal to 0. With min_count=1, it equals NaN if there are less than min_count (1) real values in the series.;
    # The difference between 0 and NaN is important for plots (NaNs) aren't shown and for statistics e.g., mean.
    return df.groupby(df.index.to_period('Y').year).sum(min_count=1)

def prepare_for_facetgrid(df):
    return df.melt(ignore_index=False).reset_index().rename(columns={'value': 'BA'})

def drop_model(df, modelname):
    return df.drop(modelname, axis=1, level=1)

def drop_models(df, modelnames):
    for modelname in modelnames:
        df = drop_model(df, modelname)
    return df

def select_model(df, modelname):
    return df.loc[:, (slice(None), modelname)]

def select_models(df, modelnames):
    return df.loc[:, (slice(None), modelnames)]

def select_region(df, regionname):
    return df.loc[:, (regionname, slice(None))]

def to_anomaly(df):
    return df - df.mean()

def to_relative_anomaly(df):
    return (df - df.mean())/df.mean()

def df_constrain_time(df, start, end):
    return df.loc[str(start):str(end)]

def df_to_global(df):
    if 'Observation' in df.columns.names:
        level = 'Observation'
    elif 'Model' in df.columns.names:
        level = 'Model'
    else:
        raise ValueError("Observation or Model should be in the column level names (df.column.names)")
    return df.groupby(level=level, axis=1).sum(min_count=1)








#######Get weights#####


obsclim_df = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/AR6_obsclim_df.pkl')
counterclim_df = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/AR6_counterclim_df.pkl')
model_weights = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/NME3_Weights.pkl')

##Set up obsclim relative anomaly
obsclim_global = obsclim_df.groupby('Model', axis=1).sum()
obsclim = df_constrain_time(obsclim_global, 2003, 2019)
model_DF = (obsclim-obsclim.mean(axis=0))/obsclim.mean(axis=0)

##Set up obs relative anomaly
    
# Make GFED5 NME
obs = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/AR6_obs_df.pkl')
obs_global = obs_df.groupby('Observation', axis=1).sum()
obs_constrainedtime = df_constrain_time(obs_global, 2003, 2019)
obs_GFED = (obs_constrainedtime['GFED5'])
obs_DF_GFED = (obs_GFED-obs_GFED.mean(axis=0))/obs_GFED.mean(axis=0)

GFED5weights = (np.abs(np.log(model_DF+1.0+(1/204))).subtract(np.log(obs_DF_GFED+1.0+(1/204)), axis=0))# Model vs Obs anomaly
GFED5weights[:] = GFED5weights.mean(axis=0)


# Make FireCCI NME
obs_CCI = (obs_constrainedtime['FireCCI5.1'])
obs_DF_CCI = (obs_CCI-obs_CCI.mean(axis=0))/obs_CCI.mean(axis=0)
CCIweights = (np.abs(np.log(model_DF+1.0+(1/204))).subtract(np.abs(np.log(obs_DF_CCI+1.0+(1/204))),axis=0))# Model vs Obs anomaly
CCIweights[:] = CCIweights.mean(axis=0)


OBSweights = (GFED5weights + CCIweights)/2
print(OBSweights)


###Add uncertainty####

import math
obsclim = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/AR6_obsclim_df.pkl')
counterclim = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/AR6_counterclim_df.pkl')

obsclim_global = obsclim_df.groupby('Model', axis=1).sum()
counterclim_global = counterclim_df.groupby('Model', axis=1).sum()
obsclim = df_constrain_time(obsclim_global, 2003, 2019)
counterclim = df_constrain_time(counterclim_global, 2003, 2019)



#######

def run_RR(obsclim, counterclim, add_noise = True):
   
    def add_random_noise(series):
        Noise = series.copy()
        for modelname in Noise.columns.unique(level='Model'):
            Noise[modelname]= (np.random.normal(loc=0, scale=np.sqrt(math.pi/2), size=len(series)))
        series = np.log(series) + Noise * OBSweights
        series = np.exp(series)
        return(series)
    
    #Calculate the uncertainty by adding noise into distribution
    if add_noise:
        obsclim = add_random_noise(obsclim) 
        counterclim = add_random_noise(counterclim) 

    #Calculate Relative anomaly and resample
    rng = np.random.default_rng(SEED)
    df = (obsclim-counterclim.mean(axis=0))/counterclim.mean(axis=0)
    series = pd.Series(df.values.ravel())
    obsclim_sampled_global = series.sample(n=100000, weights=np.tile(select_region(model_weights, 'Global').values.ravel(), len(df.index)), random_state=rng, replace=True)

    rng = np.random.default_rng(SEED)
    df = (counterclim-counterclim.mean(axis=0))/counterclim.mean(axis=0)
    series = pd.Series(df.values.ravel())
    counterclim_sampled_global = series.sample(n=100000, weights=np.tile(select_region(model_weights, 'Global').values.ravel(), len(df.index)), random_state=rng, replace=True)
    
    ##Calculate the Probability Ratio
    ALL = (np.count_nonzero(obsclim_sampled_global > (counterclim_sampled_global.mean(axis=0))))
    NAT = (np.count_nonzero(counterclim_sampled_global > (counterclim_sampled_global.mean(axis=0))))
    PR = ALL/NAT
    
    burn_change = obsclim_sampled_global.mean(axis=0)*100 - counterclim_sampled_global.mean(axis=0)*100
    return PR, burn_change


PR_noNoise, change_noNoise = run_RR(obsclim, counterclim, False)
PRuncertEnsemble = np.array([run_RR(obsclim, counterclim, True) for i in range(1000)])



#####Get the 10th and 90th percentile for the new uncertainty ensemble
PR_uncertainty_range = np.percentile(PRuncertEnsemble[:,0], [10, 90])
change_uncertainty_range = np.percentile(PRuncertEnsemble[:,1], [10, 90])

def print_rounded_result(string, result):
    print(string + str(np.round(result, 2)))
    
print_rounded_result("PR: " , PR_noNoise)
print_rounded_result("    range:", PR_uncertainty_range)
print_rounded_result("    p-value:", np.mean(PRuncertEnsemble[:,0] < 1.0))


print_rounded_result("change in burnt area:", change_noNoise)
print_rounded_result("    range:", change_uncertainty_range)
print_rounded_result("    p-value:", np.mean(PRuncertEnsemble[:,1] < 0.0))


PR: 1.43
    range:[1.36 1.42]
    p-value:0.0
change in burnt area:16.43
    range:[15.66 17.09]
    p-value:0.0


