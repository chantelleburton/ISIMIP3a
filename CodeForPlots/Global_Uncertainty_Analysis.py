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
#model_weights = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/NME3_Weights.pkl')
model_weights = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Weights_optimal_sigmaD.pkl')


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

GFED5_NME = model_DF.copy()
for modelname in model_DF.columns:
    GFED5_NME[modelname] = GFED5_NME[modelname].subtract(obs_DF['GFED5'], axis=0)
GFED5_NME = np.abs(GFED5_NME)
GFED5_NME[:] = GFED5_NME.mean(axis=0)




# Make temporal NME for FireCCI 
obs = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/AR6_obs_df.pkl')#2 obs, CCI5.1 & MCD64
obs_time = df_constrain_time(obs, 2003, 2019)
obs_global = obs_time.groupby(by='Observation', axis=1).sum(min_count=1) 
obs_CCI = obs_global[['FireCCI5.1']]
obs_DF = (obs_CCI-obs_CCI.mean(axis=0))/obs_CCI.mean(axis=0)

CCI_NME = model_DF.copy()
for modelname in model_DF.columns:
    CCI_NME[modelname] = CCI_NME[modelname].subtract(obs_DF['FireCCI5.1'], axis=0)
CCI_NME = np.abs(CCI_NME)
CCI_NME[:] = CCI_NME.mean(axis=0)


NME = (GFED5_NME + CCI_NME)/2


################################ Get Data ###################################

import math
obsclim = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/AR6_obsclim_df.pkl')
counterclim = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/AR6_counterclim_df.pkl')

obsclim_global = obsclim_df.groupby('Model', axis=1).sum()
counterclim_global = counterclim_df.groupby('Model', axis=1).sum()

model_weights = model_weights.reindex(sorted(model_weights.columns), axis = 1)
NME = NME.reindex(sorted(NME.columns), axis = 1)

def add_uncertainty(obsclim, counterclim):
   
    obsclim = obsclim.reindex(sorted(obsclim.columns), axis = 1)
    counterclim = counterclim.reindex(sorted(counterclim.columns), axis = 1)

################################ Add Uncertainty ###################################

    def run_RR(add_noise = True, i=0):
        rng = np.random.default_rng(SEED)
    
        #(2b) add the random noise
        def add_random_noise(model, NME):
            Noise = model.copy()
            for modelname in model.columns.unique(level='Model'):
                Noise[modelname] = np.random.normal(loc=0, scale=np.sqrt(math.pi/2), size=204) # Get the random numbers
            model = model + Noise * NME.to_numpy()
            return model  
    
        #(2a) add noise if TRUE
        def add_noise_to_experiment(model):
            if add_noise: model = add_random_noise(model, NME)
            return (model)

        #(4a)### find RA after adding in the noise
        def randomly_sample_RA(experiment):
            RA = (experiment-counterclim_mean)/counterclim_mean
            RA = pd.Series(RA.values.ravel())
            return np.random.choice(RA, 100000, p = weights, replace=True)    

        #(2) Call function add_noise
        obsclim_series = add_noise_to_experiment(obsclim)
        counterclim_series = add_noise_to_experiment(counterclim) 
        counterclim_mean = counterclim_series.mean(axis=0)

        #(3) get the weights ready to use
        weights = np.tile(select_region(model_weights, 'Global').values.ravel(), 204)
        weights = weights / np.sum(weights)

        #(4) make the bootstrapped, resampled relative anomaly
        obsclim_sampled_global = randomly_sample_RA(obsclim_series)
        counterclim_sampled_global = randomly_sample_RA(counterclim_series)
    
        #(5)
        ##Calculate the Probability Ratio
        ALL = (np.count_nonzero(obsclim_sampled_global > (counterclim_sampled_global.mean(axis=0))))
        NAT = (np.count_nonzero(counterclim_sampled_global > (counterclim_sampled_global.mean(axis=0))))
        PR = ALL/NAT
    
        burn_change = np.mean(obsclim_sampled_global)
    
        return PR, burn_change

    #(1) Call the functions
    PR_noNoise, change_noNoise  = run_RR(False)
    PRuncertEnsemble = np.array([run_RR(True) for i in range(1000)])

    ################################ Print results ###################################

    ##Get the 10th and 90th percentile for the new uncertainty ensemble
    PR_uncertainty_range = np.percentile(PRuncertEnsemble[:,0], [10, 90])
    change_uncertainty_range = np.percentile(PRuncertEnsemble[:,1], [10, 90])

    CHG = (np.percentile(PRuncertEnsemble[:,0], 50))
    RNG = (PR_uncertainty_range[1]) - CHG
    print('PR:',f'{CHG.round(2)} ± {RNG.round(2)}')

    CHG = np.percentile(PRuncertEnsemble[:,1], 50)*100
    RNG = (change_uncertainty_range[1]*100) - CHG
    print('BA chg:',f'{CHG.round(2)} ± {RNG.round(2)}')



###All exps PD / DHF / ALL ####
#PD
print('PD')
obsclim = df_constrain_time(obsclim_global, 2003, 2019)
counterclim = df_constrain_time(counterclim_global, 2003, 2019)
add_uncertainty(obsclim, counterclim)
#DHF
print('DHF')
obsclim = df_constrain_time(counterclim_global, 2003, 2019)
counterclim = df_constrain_time(counterclim_global, 1901, 1917)
add_uncertainty(obsclim, counterclim)
#All
print('ALL')
obsclim = df_constrain_time(obsclim_global, 2003, 2019)
counterclim = df_constrain_time(counterclim_global, 1901, 1917)
add_uncertainty(obsclim, counterclim)





'''
#Results

PD
1.5 ± 0.06
21.97 ± 0.2
DHF
0.34 ± 0.03
-22.18 ± 0.11
ALL
0.78 ± 0.04
-6.32 ± 0.13



'''






