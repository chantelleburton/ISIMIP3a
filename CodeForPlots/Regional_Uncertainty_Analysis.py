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
def select_model(df, modelname):
    return df.loc[:, (slice(None), modelname)]
SEED = int.from_bytes("They're taking the hobbits to Isengard!".encode('utf-8'), byteorder='big')




#######################################Get weights########################################################

obsclim_df = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/AR6_obsclim_df.pkl')
counterclim_df = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/AR6_counterclim_df.pkl')
model_weights = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/NME3_Weights.pkl')

##Set up obsclim relative anomaly
obsclim = df_constrain_time(obsclim_df, 2003, 2019)
model_DF = (obsclim-obsclim.mean(axis=0))/obsclim.mean(axis=0)

##Set up obs relative anomaly

    
# Make GFED5 NME
obs = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/AR6_obs_df.pkl')
obs_GFED = df_constrain_time(select_model(obs, 'GFED5'), 2003, 2019)
obs_DF_GFED = (obs_GFED-obs_GFED.mean(axis=0))/obs_GFED.mean(axis=0)
#print(obs_DF_GFED)

GFED5weights = (np.log(model_DF+1.0+(1/204))).subtract(np.log(obs_DF_GFED+1.0+(1/204))).droplevel('Observation', axis=1)# Model vs Obs anomaly
GFED5weights = np.abs(GFED5weights)
GFED5weights[:] = GFED5weights.mean(axis=0)
#print(GFED5weights)

# Make FireCCI NME
obs = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/AR6_obs_df.pkl')
obs_CCI = df_constrain_time(select_model(obs, 'FireCCI5.1'), 2003, 2019)
obs_DF_CCI = (obs_CCI-obs_CCI.mean(axis=0))/obs_CCI.mean(axis=0)

FireCCIweights = (np.log(model_DF+1.0+(1/204))).subtract(np.log(obs_DF_CCI+1.0+(1/204))).droplevel('Observation', axis=1)# Model vs Obs anomaly
FireCCIweights = (np.abs(FireCCIweights))
FireCCIweights[:] = FireCCIweights.mean(axis=0)
#print(FireCCIweights)

#Average both obs datasets
OBSweights = (GFED5weights + FireCCIweights)/2
#print(OBSweights)


####################################################### Load data######################################################
import math
obsclim = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/AR6_obsclim_df.pkl')
counterclim = pd.read_pickle(f'/scratch/cburton/scratch/ISIMIP3a/Data/AR6_counterclim_df.pkl')

obsclim = df_constrain_time(obsclim, 2003, 2019)
counterclim = df_constrain_time(counterclim, 2003, 2019)

obsclim = obsclim.reindex(sorted(obsclim.columns), axis = 1)
counterclim = counterclim.reindex(sorted(counterclim.columns), axis = 1)
model_weights = model_weights.reindex(sorted(model_weights.columns), axis = 1)
OBSweights = OBSweights.reindex(sorted(OBSweights.columns), axis = 1)
######################################################## Get results ###################################################

### Regional Uncertainty results

#######
def run_RR(obsclim, counterclim, add_noise = True, i = 0):
    
    def add_random_noise(model):
        Noise = model.copy()
        for regionname in model.columns.unique(level='Region'):                 
            for modelname in model.columns.unique(level='Model'):
                 Noise[regionname,modelname]= (np.random.normal(loc=0, scale=np.sqrt(math.pi/2), size=len(model)))
        log_model = np.log(model+(1/len(model)))
        series = log_model + Noise * OBSweights
        series = series - np.mean(series, axis=0)
        series = series / np.std(series, axis=0)
        series = series * np.std(log_model)/np.std(series)
        series = series + np.mean(log_model)
        series = np.exp(series)
        return(series)
    
    #Add noise second if TRUE (calculate the uncertainty by adding noise into distribution)
    if add_noise:
        obsclim = add_random_noise(obsclim) 
        counterclim = add_random_noise(counterclim) 

    #Calculate Relative anomaly and resample
    rng = np.random.default_rng(SEED)
    df = (obsclim-counterclim.mean(axis=0))/counterclim.mean(axis=0)
    obsclim_sampled_regional = pd.DataFrame(index=np.arange(100000), columns=df.columns.unique(level='Region'))
    for regionname in df.columns.unique(level='Region'):
        region_series = pd.Series(df.loc[slice(None), regionname].values.ravel()) 
        obsclim_sampled_regional.loc[slice(None), regionname] = region_series.sample(n=100000, weights=np.tile(select_region(model_weights, regionname).values.ravel(), len(df.index)), random_state=rng, replace=True).values
    obsclim_sampled_regional.replace([np.inf, -np.inf], np.nan, inplace=True)
    obsclim_sampled_regional.dropna(inplace=True)

    rng = np.random.default_rng(SEED)
    df = (counterclim-counterclim.mean(axis=0))/counterclim.mean(axis=0)
    counterclim_sampled_regional = pd.DataFrame(index=np.arange(100000), columns=df.columns.unique(level='Region'))
    for regionname in df.columns.unique(level='Region'):
        region_series = pd.Series(df.loc[slice(None), regionname].values.ravel())
        counterclim_sampled_regional.loc[slice(None), regionname] = region_series.sample(n=100000, weights=np.tile(select_region(model_weights, regionname).values.ravel(), len(df.index)), random_state=rng, replace=True).values
    counterclim_sampled_regional.replace([np.inf, -np.inf], np.nan, inplace=True)
    counterclim_sampled_regional.dropna(inplace=True)


    ##Calculate the Probability Ratio
    ALL = (np.count_nonzero(obsclim_sampled_regional > (counterclim_sampled_regional.mean(axis=0))))
    NAT = (np.count_nonzero(counterclim_sampled_regional > (counterclim_sampled_regional.mean(axis=0))))
    PR = ALL/NAT

    burn_change = obsclim_sampled_regional.mean(axis=0)*100 - counterclim_sampled_regional.mean(axis=0)*100


    PRlist = []
    for regionname in df.columns.unique(level='Region'):
        ALL = (np.count_nonzero(obsclim_sampled_regional.loc[slice(None), regionname] > (counterclim_sampled_regional.loc[slice(None), regionname].mean(axis=0))))
        NAT = (np.count_nonzero(counterclim_sampled_regional.loc[slice(None), regionname] > (counterclim_sampled_regional.loc[slice(None), regionname].mean(axis=0))))
        PR = ALL/NAT
        PRlist.append(PR)
    
    if add_noise:
        return burn_change, PRlist

    
    if add_noise == False:
        Table = pd.DataFrame(burn_change, columns=['BurnChange'])
        Table['PR'] = PRlist
        return Table


### FIRST read in the data###
Table = (run_RR(obsclim, counterclim, False))

Burn, PR = zip(*[run_RR(obsclim, counterclim, True) for i in range(1000)])

PRUncertTable = pd.DataFrame(PR, columns=obsclim.columns.unique(level='Region'))
BurnUncertTable = pd.DataFrame(Burn, columns=obsclim.columns.unique(level='Region'))

Table["Burn10"] = (BurnUncertTable[:].quantile(0.1))
Table["Burn90"] = (BurnUncertTable[:].quantile(0.9))
Table["PR10"] = (PRUncertTable[:].quantile(0.1))
Table["PR90"] = (PRUncertTable[:].quantile(0.9))
Table = Table.round(decimals=2)

Table = Table[['BurnChange', 'Burn10', 'Burn90', 'PR', 'PR10', 'PR90']]
print(Table.astype(float).round(2))

Table.to_pickle("~/GitHub/ISIMIP3a/Supplementary_Data.Uncertainty2.pkl") 
print('Saved')
