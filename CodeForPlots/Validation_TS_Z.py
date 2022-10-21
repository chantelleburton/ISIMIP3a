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

def CollapseToTimeseries(cube):
    coords = ('longitude', 'latitude')
    for coord in coords:
        if not cube.coord(coord).has_bounds():
            cube.coord(coord).guess_bounds()
    weights = iris.analysis.cartography.area_weights(cube)
    cube = cube.collapsed(coords, iris.analysis.SUM, weights = weights)/1E12*100/100 # Convert to Mkm2, and then to Mha, and change from percent to frac 
    return cube

def MakeAnomaly(cube):
    Mean = cube.collapsed('time', iris.analysis.MEAN)
    Anomaly = cube-Mean
    StdDev = np.std(cube.data)
    StdDev = Anomaly/StdDev
    return StdDev

folder = '/scratch/cburton/scratch/ISIMIP3a/Data/'
d = {}

##### Load Model Data ##### 
models = ['jules', 'classic','visit', 'ssib4', 'LPJ-GUESS-SPITFIRE', 'LPJ-GUESS-SIMFIRE']
for model in models:
    print(model)
    cube = iris.load_cube(folder+model+'*_gswp3-w5e5_obsclim_histsoc_default_burntarea-total_global_monthly_1901_2019.nc')
    if model == 'LPJ-GUESS-SPITFIRE' or model == 'LPJ-GUESS-SIMFIRE':
        cube.data = ma.masked_where(np.isnan(cube.data),cube.data)
    else:
       cube = UpdateTime(cube)
    cube = MonthlyToAnnual(cube)
    cube = ConstrainTime(cube)
    if model == 'classic':
        cube = cube*100 #Convert frac to percent
    if model == 'ssib4':
        cube = cube*30 #Convert %/day to %/month
    cube = CollapseToTimeseries(cube) 
    d[(model)] = MakeAnomaly(cube)

##### Load Observations #####  
obs = ['GFED4.1s', 'FireCCI5.1', 'FireCCILT11']
for ob in obs:
    print (ob)
    cube = iris.load_cube(folder+ob+'_Burned_Percentage.nc')
    cube = MonthlyToAnnual(cube)
    cube = ConstrainTime(cube)
    cube = CollapseToTimeseries(cube)
    d[(ob)] = MakeAnomaly(cube)


##### Make Timeseries Plot ######
AllData = ['jules', 'classic','visit', 'ssib4', 'LPJ-GUESS-SPITFIRE', 'LPJ-GUESS-SIMFIRE','GFED4.1s', 'FireCCI5.1', 'FireCCILT11']
colors= ('white','lightblue', 'coral', 'sandybrown', 'lightgreen', 'plum', 'pink', 'black', 'blue', 'darkblue') 
n = 1
x = np.arange(2001,2017)
for Data in AllData: 
    slopen, interceptn, r_valuen, p_valuen, std_errn = stats.linregress(x,d[Data].data)
    plt.plot(x,d[Data].data, label=Data, color=colors[n % len(colors)])
    plt.plot(x, interceptn + slopen*x, color=colors[n % len(colors)])
    n = n+1

plt.ylabel('Z-score')
plt.title('Annual Total Burnt Area Z-score')
plt.legend(ncol=2, loc='best')
plt.show()


##### Make Range Plot ######
Models = (d['classic'].data,d['jules'].data,d['visit'].data,d['ssib4'].data,d['LPJ-GUESS-SPITFIRE'].data)
ModelMax=np.amax(Models, axis=0)
ModelMin=np.amin(Models, axis=0)

Obs = (d['GFED4.1s'].data,d['FireCCI5.1'].data,d['FireCCILT11'].data)
ObsMax=np.amax(Obs, axis=0)
ObsMin=np.amin(Obs, axis=0)

plt.fill_between(x, ModelMin, ModelMax, color='lightgrey', label='Model range')
plt.fill_between(x, ObsMin, ObsMax, color='grey', label='Obs range')
plt.legend()
plt.ylabel('Z-score')
plt.title('Annual total burned area Z-score')
plt.show()



