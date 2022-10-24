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
import statsmodels.api as sm
from scipy.stats import norm, uniform
import cartopy.io.shapereader as shpreader
import ascend
from ascend import shape



def MaskRegion(cube, region):
    regions = ascend.shape.load_shp(ascend.EXAMPLE_GIORGI, short_name=region)
    domain_shape = regions.unary_union()   
    MaskedCube = domain_shape.mask_cube(cube)
    return MaskedCube

def UpdateTime(cube):
    timeco = cube.coord('time')
    assert timeco.units == cf_units.Unit('months since 1901-01-01', calendar='360_day')
    timeco.units = cf_units.Unit('days since 1901-01-01', calendar='360_day')
    timeco.points = timeco.points * 30.
    return cube

def ConstrainTime(cube, date1, date2):
    date = iris.Constraint(time=lambda cell: date1 <= cell.point.year <= date2)
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


##Change region here
region = 'SSA'

##### Load Model Data ##### 
folder = '/scratch/cburton/scratch/ISIMIP3a/Data/'
models = ['jules', 'classic','visit', 'ssib4', 'LPJ-GUESS-SPITFIRE', 'LPJ-GUESS-SIMFIRE']
obs = ['GFED', 'CCI51', 'CCILT11']
d = {}
for model in models:
    for ob in obs:
        cube = iris.load_cube(folder+model+'*_gswp3-w5e5_obsclim_histsoc_default_burntarea-total_global_monthly_1901_2019.nc')
        if model == 'LPJ-GUESS-SPITFIRE' or model == 'LPJ-GUESS-SIMFIRE':
            cube.data = ma.masked_where(np.isnan(cube.data),cube.data)
        else:
            cube = UpdateTime(cube)
        cube = MonthlyToAnnual(cube)
        if model == 'classic':
            cube = cube*100 #Convert frac to percent
        if model == 'ssib4':
            cube = cube*30 #Convert %/day to %/month
        if ob == 'GFED':
            date1=1997
            date2=2016
        elif ob == 'CCI51':
            date1=2001
            date2=2019
        elif ob == 'CCILT11':
            date1=1982
            date2=2017
        cube = ConstrainTime(cube, date1, date2)
        #cube = MaskRegion(cube, region=region)
        cube = CollapseToTimeseries(cube)
        d[(model+ob)] = MakeAnomaly(cube)


##### Load Observations #####  
obs = ['GFED4.1s', 'FireCCI5.1', 'FireCCILT11']
for ob in obs:
    cube = iris.load_cube(folder+ob+'_Burned_Percentage.nc')
    cube = MonthlyToAnnual(cube)
    if ob == 'GFED4.1s':
        date1=1997
        date2=2016
    elif ob == 'FireCCI5.1':
        date1=2001
        date2=2019
    elif ob == 'FireCCILT11':
        date1=1982
        date2=2018
    cube = ConstrainTime(cube, date1, date2)
    #cube = MaskRegion(cube, region=region)
    cube = CollapseToTimeseries(cube)
    d[(ob)] = MakeAnomaly(cube)


# Quantile-quantile / Scatter plot
plt.figure()
plt.plot(np.sort(d['GFED4.1s'].data), np.sort(d['julesGFED'].data), marker='x', color='lightblue', label='JULES v GFED')
plt.plot(np.sort(d['GFED4.1s'].data), np.sort(d['classicGFED'].data), marker='x', color='coral', label='CLASSIC v GFED')
plt.plot(np.sort(d['GFED4.1s'].data), np.sort(d['visitGFED'].data), marker='x', color='sandybrown', label='VISIT v GFED')
plt.plot(np.sort(d['GFED4.1s'].data), np.sort(d['ssib4GFED'].data), marker='x', color='lightgreen', label='SSIB4 v GFED')
plt.plot(np.sort(d['GFED4.1s'].data), np.sort(d['LPJ-GUESS-SPITFIREGFED'].data), marker='x', color='plum', label='SPITFIRE v GFED')

plt.plot(np.sort(d['FireCCI5.1'].data), np.sort(d['julesCCI51'].data), color='lightblue', label='CCI51')
plt.plot(np.sort(d['FireCCI5.1'].data), np.sort(d['classicCCI51'].data), color='coral')
plt.plot(np.sort(d['FireCCI5.1'].data), np.sort(d['visitCCI51'].data), color='sandybrown')
plt.plot(np.sort(d['FireCCI5.1'].data), np.sort(d['ssib4CCI51'].data), color='lightgreen')
plt.plot(np.sort(d['FireCCI5.1'].data), np.sort(d['LPJ-GUESS-SPITFIRECCI51'].data), color='plum')

plt.plot(np.sort(d['FireCCILT11'].data), np.sort(d['julesCCILT11'].data), linestyle='dashed', color='lightblue', label='CCI11')
plt.plot(np.sort(d['FireCCILT11'].data), np.sort(d['classicCCILT11'].data), linestyle='dashed', color='coral')
plt.plot(np.sort(d['FireCCILT11'].data), np.sort(d['visitCCILT11'].data), linestyle='dashed', color='sandybrown')
plt.plot(np.sort(d['FireCCILT11'].data), np.sort(d['ssib4CCILT11'].data), linestyle='dashed', color='lightgreen')
plt.plot(np.sort(d['FireCCILT11'].data), np.sort(d['LPJ-GUESS-SPITFIRECCILT11'].data), linestyle='dashed', color='plum')



#For absolute value plot
#axis_min=350
#axix_max=850
#For std dev plot
axis_min=-2.5
axix_max=3.5
#For anomaly plot
#axis_min=-120
#axix_max=200

plt.xlabel('Observed burned area Z-score')
plt.ylabel('FireMIP burned area Z-score')
#plt.xlabel('Observed burned area anomaly (Mha)')
#plt.ylabel('FireMIP burned area anomaly (Mha)')
plt.plot([axis_min,axix_max],[axis_min,axix_max],color="k")
plt.ylim(axis_min,axix_max)
plt.xlim(axis_min,axix_max)
plt.legend()
plt.show()







'''
TEST : using qqplot is exactly the same as above

import statsmodels.api as sm
from statsmodels.graphics.gofplots import qqplot_2samples
x = GFED.data
y = d['classicGFED'].data
pp_x = sm.ProbPlot(x)
pp_y = sm.ProbPlot(y)
qqplot_2samples(pp_x, pp_y)
plt.scatter(np.sort(GFED.data), np.sort(d['classicGFED'].data), marker='x', color='red', label='JULES v GFED')
axis_min=-80
axix_max=60
plt.plot([axis_min,axix_max],[axis_min,axix_max],color="k")
plt.ylim(axis_min,axix_max)
plt.xlim(axis_min,axix_max)
plt.xlabel('Observed global total burned area (Mha)')
plt.ylabel('FireMIP global total burned area (Mha)')
plt.show()
exit()

'''


'''
#Alternative Q-Q plot, opens two separate plots

from statsmodels.graphics.gofplots import qqplot_2samples
#https://www.statsmodels.org/stable/generated/statsmodels.graphics.gofplots.qqplot_2samples.html#statsmodels.graphics.gofplots.qqplot_2samples
x = FireCCI51.data

y1 = d['jules'].data
pp_x = sm.ProbPlot(x)
pp_y1 = sm.ProbPlot(y1)
qqplot_2samples(pp_x, pp_y1)

y2 = d['classic'].data
pp_y2 = sm.ProbPlot(y2)
qqplot_2samples(pp_x, pp_y2)

plt.xlabel('FireCCI51 global total burned area (Mha)')
plt.ylabel('FireMIP global total burned area (Mha)')
plt.plot([min(x),max(x)],[min(x),max(x)],color="k")
plt.show()
raise Exception('stop')
'''








