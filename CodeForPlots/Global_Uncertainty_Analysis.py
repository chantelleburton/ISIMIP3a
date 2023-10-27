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

    
def log_diff(obs_DF, model_DF):
    return np.mean(np.abs(np.log(model_DF+1.0) - np.log(obs_DF[:,0]+1.0)))

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

GFED5weights = []
for modelname in model_DF.columns:    
    Weights = log_diff(obs_DF.values, model_DF[modelname].values)
    GFED5weights.append(Weights)
print(GFED5weights)

# Make temporal NME for FireCCI 
obs = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/AR6_obs_df.pkl')#2 obs, CCI5.1 & MCD64
obs_time = df_constrain_time(obs, 2003, 2019)
obs_global = obs_time.groupby(by='Observation', axis=1).sum(min_count=1) 
obs_CCI = obs_global[['FireCCI5.1']]
obs_DF = (obs_CCI-obs_CCI.mean(axis=0))/obs_CCI.mean(axis=0)

FireCCIweights = []
for modelname in model_DF.columns:    
    Weights = log_diff(obs_DF.values, model_DF[modelname].values)
    FireCCIweights.append(Weights)
print(FireCCIweights)

#"Broadcast" the results for GFED into a big dataframe to get the right size later (204 rows)
GFED5weights_df = pd.DataFrame(columns=obsclim.columns.unique(level='Model'))
for n in np.arange(204):
    GFED5weights_df.loc[len(GFED5weights_df)] = GFED5weights

OBSweights_df = (GFED5weights_df + FireCCIweights_df) / 2
print(OBSweights_df)




################################ Get Data ###################################

import math
obsclim = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/AR6_obsclim_df.pkl')
counterclim = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/AR6_counterclim_df.pkl')

obsclim_global = obsclim_df.groupby('Model', axis=1).sum()
counterclim_global = counterclim_df.groupby('Model', axis=1).sum()
obsclim = df_constrain_time(obsclim_global, 2003, 2019)
counterclim = df_constrain_time(counterclim_global, 2003, 2019)


################################ Add Uncertainty ###################################

def run_RR(add_noise = True):
    rng = np.random.default_rng(SEED)
    
    def add_random_noise(model, noise):
    
        # Your logic here
        random_numbers = np.random.normal(loc=0, scale=np.sqrt(math.pi/2), size=len(model)) # Get the random numbers
        log_model = np.log(model)
        series = log_model + random_numbers * noise.values
        series = series - np.mean(series)
        series = series * np.std(log_model)/np.std(series)
        series = series + np.mean(log_model)
        series = np.exp(series)
        return series  # Example: Adding the two columns
    
    def add_noise_to_experiment(models):
        if add_noise: models = models.apply(lambda col1: add_random_noise(col1, OBSweights_df[col1.name]), axis=0)
        
        return pd.Series(models.values.ravel())
    
    #Calculate the uncertainty for obsclim
    obsclim_series = add_noise_to_experiment(obsclim)
    counterclim_series = add_noise_to_experiment(counterclim) 
    counterclim_mean = counterclim_series.mean(axis=0)
    
    weights = np.tile(select_region(model_weights, 'Global').values.ravel(), 204)
    weights = weights / np.sum(weights)
    
    def randomly_sample_RA(experiment):
        RA = (experiment-counterclim_mean)/counterclim_mean
        return np.random.choice(RA, 10000, p = weights, replace=True)
    
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
    range:[1.4  1.48]
    p-value:0.0
change on burnt area:0.16
    range:[0.15 0.16]
    p-value:0.0

'''


