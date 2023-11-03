import os
import numpy as np
import seaborn as sns
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
warnings.filterwarnings('ignore')


def df_constrain_time(df, start, end):
    return df.loc[str(start):str(end)]
def select_region(df, regionname):
    return df.loc[:, (regionname, slice(None))]
SEED = int.from_bytes("They're taking the hobbits to Isengard!".encode('utf-8'), byteorder='big')


obsclim_df = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/AR6_obsclim_df.pkl')
counterclim_df = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/AR6_counterclim_df.pkl')
model_weights = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/NME3_Weights.pkl')



#######################################Get weights########################################################


#Set up obsclim relative anomaly
obsclim_global = obsclim_df.groupby('Model', axis=1).sum()
obsclim = df_constrain_time(obsclim_global, 2003, 2019)
model_DF = (obsclim-obsclim.mean(axis=0))/obsclim.mean(axis=0)

#Set up obs relative anomaly

# Make temporal NME for GFED5 first
obs = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/AR6_obs_df.pkl')#2 obs, CCI5.1 & MCD64
obs_time = df_constrain_time(obs, 2003, 2019)
obs_global = obs_time.groupby(by='Observation', axis=1).sum(min_count=1) 
obs_GFED = obs_global[['GFED5']]
obs_DF = (obs_GFED-obs_GFED.mean(axis=0))/obs_GFED.mean(axis=0)


GFED5weights = model_DF.copy()
for modelname in model_DF.columns:
    GFED5weights[modelname] = np.log(GFED5weights[modelname]+1.0+(1/204)).subtract(np.log(obs_DF['GFED5']+1.0+(1/204)), axis=0)
GFED5weights = np.abs(GFED5weights)
GFED5weights[:] = GFED5weights.mean(axis=0)


# Make temporal NME for FireCCI 
obs = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/AR6_obs_df.pkl')#2 obs, CCI5.1 & MCD64
obs_time = df_constrain_time(obs, 2003, 2019)
obs_global = obs_time.groupby(by='Observation', axis=1).sum(min_count=1) 
obs_CCI = obs_global[['FireCCI5.1']]
obs_DF = (obs_CCI-obs_CCI.mean(axis=0))/obs_CCI.mean(axis=0)


CCIweights = model_DF.copy()
for modelname in model_DF.columns:
    CCIweights[modelname] = np.log(CCIweights[modelname]+1.0+(1/204)).subtract(np.log(obs_DF['FireCCI5.1']+1.0+(1/204)), axis=0)
CCIweights = np.abs(CCIweights)
CCIweights[:] = CCIweights.mean(axis=0)


OBSweights = (GFED5weights + CCIweights)/2


################################ Get Data ###################################

import math
obsclim = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/AR6_obsclim_df.pkl')
counterclim = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/AR6_counterclim_df.pkl')

obsclim_global = obsclim_df.groupby('Model', axis=1).sum()
counterclim_global = counterclim_df.groupby('Model', axis=1).sum()
obsclim = df_constrain_time(obsclim_global, 2003, 2019)
counterclim = df_constrain_time(counterclim_global, 2003, 2019)

obsclim = obsclim.reindex(sorted(obsclim.columns), axis = 1)
counterclim = counterclim.reindex(sorted(counterclim.columns), axis = 1)
model_weights = model_weights.reindex(sorted(model_weights.columns), axis = 1)
OBSweights = OBSweights.reindex(sorted(OBSweights.columns), axis = 1)

################################ Add Uncertainty ###################################

def run_RR(add_noise = True, i=0):
    rng = np.random.default_rng(SEED)
    
    def add_random_noise(model, OBSweights):
        Noise = model.copy()
        for modelname in model.columns.unique(level='Model'):
            Noise[modelname] = np.random.normal(loc=0, scale=np.sqrt(math.pi/2), size=len(model)) # Get the random numbers
        log_model = np.log(model+(1/len(model)))
        series = log_model + Noise * OBSweights
        series = series - np.mean(series, axis=0)
        series = series / np.std(series, axis=0)
        series = series * np.std(log_model)/np.std(series)
        series = series + np.mean(log_model)
        series = np.exp(series)
        return series  
    
    def add_noise_to_experiment(models):
        if add_noise: models = add_random_noise(models, OBSweights)
        return (models)
    
    #Calculate the uncertainty for obsclim
    obsclim_series = add_noise_to_experiment(obsclim)
    counterclim_series = add_noise_to_experiment(counterclim) 
    counterclim_mean = counterclim_series.mean(axis=0)
    
    weights = np.tile(select_region(model_weights, 'Global').values.ravel(), 204)
    weights = weights / np.sum(weights)
    
    def randomly_sample_RA(experiment):
        RA = (experiment-counterclim_mean)/counterclim_mean
        RA = pd.Series(RA.values.ravel())
        return np.random.choice(RA, 100000, p = weights, replace=True)
    
    obsclim_sampled_global = randomly_sample_RA(obsclim_series)
    counterclim_sampled_global = randomly_sample_RA(counterclim_series)
    
    ##Calculate the Probability Ratio
    ALL = (np.count_nonzero(obsclim_sampled_global > (counterclim_sampled_global.mean(axis=0))))
    NAT = (np.count_nonzero(counterclim_sampled_global > (counterclim_sampled_global.mean(axis=0))))
    PR = ALL/NAT
    
    burn_change = np.mean(obsclim_sampled_global)
    
    return PR, burn_change


PR_noNoise, change_noNoise  = run_RR(False)
PRuncertEnsemble = np.array([run_RR(True) for i in range(1000)])

################################ Print results ###################################

##Get the 10th and 90th percentile for the new uncertainty ensemble
PR_uncertainty_range = np.percentile(PRuncertEnsemble[:,0], [10, 90])
change_uncertainty_range = np.percentile(PRuncertEnsemble[:,1], [10, 90])

def print_rounded_result(string, result):
    print(string + str(np.round(result, 2)))
    
print_rounded_result("PR: " , PR_noNoise)
print_rounded_result("    range:", PR_uncertainty_range)
print_rounded_result("    p-value:", np.mean(PRuncertEnsemble[:,0] < 1.0))


print_rounded_result("change on burnt area:", change_noNoise)
print_rounded_result("    range:", change_uncertainty_range)
print_rounded_result("    p-value:", np.mean(PRuncertEnsemble[:,1] < 0.0))



'''
PR: 1.43
    range:[1.39  1.47]
    p-value:0.0
change on burnt area:0.16
    range:[0.16 0.17]
    p-value:0.0

'''


