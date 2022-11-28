import iris
import matplotlib as mpl
import numpy as np
from iris.plot import pcolormesh
import iris.analysis.cartography
import iris.plot as iplt
from iris.coord_systems import GeogCS
import cartopy.crs as ccrs
import cf_units
import cftime
import iris.coord_categorisation
import matplotlib
import matplotlib.pyplot as plt
import iris.quickplot as qplt
from scipy import stats
import numpy.ma as ma
import matplotlib.colors as mcols



def UpdateTime(cube):
    timeco = cube.coord('time')
    assert timeco.units == cf_units.Unit('months since 1901-01-01', calendar='360_day')
    timeco.units = cf_units.Unit('days since 1901-01-01', calendar='360_day')
    timeco.points = timeco.points * 30.
    return cube

def ConstrainTime(cube):
    date = iris.Constraint(time=lambda cell: 2001 <= cell.point.year <= 2016)
    cube = cube.extract(date)
    return cube

def MonthlyToAnnual(cube):
    iris.coord_categorisation.add_year(cube, 'time', name='year')
    cube = cube.aggregated_by(['year'],iris.analysis.SUM)
    return cube

def MakeTrendMap(cube):
    TrendCube = cube.collapsed('time', iris.analysis.MEAN)
    #x = np.arange(2001,2017)
    x = np.arange(0,192)
    for i in range(len(cube.coord('longitude').points)):
        for j in range(len(cube.coord('latitude').points)):
            y = cube.data[:, j, i]
            slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
            TrendCube.data[j, i] = slope
    return TrendCube

'''
################### MAKE TREND CUBES #########################
##### Load Model Data ##### 
folder = '/scratch/cburton/scratch/ISIMIP3a/Data/'
models = ['jules', 'classic', 'visit', 'ssib4', 'LPJ-GUESS-SPITFIRE']
d = {}
for model in models:
    cube = iris.load_cube(folder+model+'*_gswp3-w5e5_obsclim_histsoc_default_burntarea-total_global_monthly_1901_2019.nc')
    if model == 'LPJ-GUESS-SPITFIRE' or model == 'LPJ-GUESS-SIMFIRE':
        cube.data = ma.masked_where(np.isnan(cube.data),cube.data)
    else:
       cube = UpdateTime(cube)
    #cube = MonthlyToAnnual(cube)
    cube = ConstrainTime(cube)
    if model == 'classic':
        cube = cube*100 #Convert frac to percent
    if model == 'ssib4':
        cube = cube*30 #Convert %/day to %/month
    d[(model)] = MakeTrendMap(cube)
    iris.save(d[(model)],'/scratch/cburton/scratch/ISIMIP3a/October/'+model+'_TrendCube_Monthly.nc') 

##### Load Observations #####   
obs = ['GFED4.1s', 'FireCCI5.1', 'FireCCILT11']
for ob in obs:
    print (ob)
    cube = iris.load_cube(folder+ob+'_Burned_Percentage.nc')
    #cube = MonthlyToAnnual(cube)
    cube = ConstrainTime(cube)
    d[(ob)] = MakeTrendMap(cube)
    iris.save(d[(ob)],'/scratch/cburton/scratch/ISIMIP3a/October/'+ob+'_TrendCube_Monthly.nc')
################### MAKE TREND CUBES #########################
'''


################### Save time by loading in pre-saved trend cubes ###################
## Make Plots ###
#minimum_log_level = 0.02
#maximum_scale_level = 2.0

minimum_log_level = 0.0001
maximum_scale_level = 0.01

cmap='RdBu_r'
anom_norm = mcols.SymLogNorm(
    linthresh=minimum_log_level,
    linscale=0.01,
    vmin=-maximum_scale_level,
    vmax=maximum_scale_level)

#tick_levels = [-2, -0.2, 0.0, 0.2, 2]
tick_levels = [-0.02, -0.002, 0.0, 0.002, 0.02]


AllData = ['GFED4.1s', 'FireCCI5.1', 'FireCCILT11','jules', 'classic','visit', 'ssib4', 'LPJ-GUESS-SPITFIRE']
n=1
for Data in AllData:
    cube = iris.load_cube('/scratch/cburton/scratch/ISIMIP3a/October/'+Data+'_TrendCube_Monthly.nc')
    plt.subplot(2,4,n)
    mesh = iplt.pcolormesh(cube, cmap=cmap, norm=anom_norm)
    bar = plt.colorbar(mesh, ticks=tick_levels, orientation="horizontal")
    bar.set_ticklabels(tick_levels)
    plt.gca().coastlines()
    plt.gcf().set_size_inches(15, 8)
    plt.title(Data)
    if n == 6:
        bar.set_label('Log scale trend in burnt area fraction (% yr $^{-2}$)')
    n = n+1
plt.show()




