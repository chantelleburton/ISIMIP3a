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
import seaborn as sns
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

def ConstrainTime(cube):
    date = iris.Constraint(time=lambda cell: 2001 <= cell.point.year <= 2016)
    cube = cube.extract(date)
    return cube
    
def MonthlyToAnnual(cube):
    iris.coord_categorisation.add_year(cube, 'time', name='year')
    cube = cube.aggregated_by(['year'],iris.analysis.SUM)
    return cube

def CollapseBySpace(cube):
    coords = ('longitude', 'latitude')
    for coord in coords:
        if not cube.coord(coord).has_bounds():
            cube.coord(coord).guess_bounds()
    weights = iris.analysis.cartography.area_weights(cube)
    cube = cube.collapsed(coords, iris.analysis.SUM, weights = weights)/1E12*100/100 # Convert to Mkm2 and then to Mha  
    return cube

def MakeAnomaly(cube):
    Mean = cube.collapsed('time', iris.analysis.MEAN)
    Anomaly = cube-Mean
    StdDev = np.std(Anomaly.data)
    StdDev = Anomaly/StdDev
    return StdDev

def CollapseByTime(cube):
    cube = cube.collapsed('time', iris.analysis.MEAN) 
    return cube

##Change region here
region = 'SSA'

##### Load Data ##### 
folder = '/scratch/cburton/scratch/ISIMIP3a/Data/'
AllData = ['GFED4.1s', 'FireCCI5.1', 'FireCCILT11','jules', 'classic','visit', 'ssib4', 'LPJ-GUESS-SPITFIRE']
colors= ('white','black', 'blue', 'darkblue','lightblue', 'coral', 'sandybrown', 'lightgreen', 'plum') 
d = {}
n = 1
for Data in AllData:
    cube = iris.load_cube(folder+Data+'*.nc')
    if Data == 'LPJ-GUESS-SPITFIRE' or Data == 'LPJ-GUESS-SIMFIRE':
        cube.data = ma.masked_where(np.isnan(cube.data),cube.data)
    elif Data == 'GFED4.1s' or Data == 'FireCCI5.1' or Data == 'FireCCILT11':
        cube = cube
    else:
        cube = UpdateTime(cube)
    cube = MonthlyToAnnual(cube)
    cube = ConstrainTime(cube)
    if Data == 'classic':
        cube = cube*100 #Convert frac to percent
    if Data == 'ssib4':
        cube = cube*30 #Convert %/day to %/month
    #cube = MaskRegion(cube, region=region)
    cube = CollapseBySpace(cube)
    d[(Data)] = MakeAnomaly(cube)
    sns.distplot(d[Data].data, hist=False, kde=True, 
                 bins=int(180/5), color = colors[n % len(colors)], 
                 hist_kws={'edgecolor':colors[n % len(colors)]},
                 kde_kws={'linewidth': 4, 'color': colors[n % len(colors)], 'label':Data})
    n = n+1
plt.legend(ncol=1, loc='best')
plt.xlabel('Normalised Burnt Area')
plt.ylabel('Probability Density')
plt.show()







'''
##Plot style 1
sns.distplot(GFED.data, kde=False,fit=stats.genextreme, color="grey",label='GFED',fit_kws={"linewidth":2.5,"color":"black"})
sns.distplot(FireCCI51.data, kde=False,fit=stats.genextreme, color="blue",label='CCI51',fit_kws={"linewidth":2.5,"color":"darkblue"})
sns.distplot(FireCCILT11.data, kde=False,fit=stats.genextreme, color="darkblue",label='CCILT11-NAT',fit_kws={"linewidth":2.5,"color":"black"})
sns.distplot(d['jules'].data, kde=False,fit=stats.genextreme, color="red",label='JULES',fit_kws={"linewidth":2.5,"color":"darkred"})
sns.distplot(d['classic'].data, kde=False,fit=stats.genextreme, color="orange",label='CLASSIC',fit_kws={"linewidth":2.5,"color":"darkorange"})
sns.distplot(d['visit'].data, kde=False,fit=stats.genextreme, color="darkred",label='VISIT',fit_kws={"linewidth":2.5,"color":"black"})
sns.distplot(d['ssib4'].data, kde=False,fit=stats.genextreme, color="yellowgreen",label='SSIB4',fit_kws={"linewidth":2.5,"color":"green"})

plt.legend()
plt.xlabel('Burnt Area (Mha)')
plt.ylabel('Probability Density')
plt.show()


##Plot style 2
plt.hist(GFED.data, color='black', alpha=0.5, label='GFED')
plt.hist(FireCCI51.data, color='blue', alpha=0.5, label='CCI51')
plt.hist(FireCCILT11.data, color='darkblue', alpha=0.5, label='CCILT11')
plt.hist(d['jules'].data, color='red', alpha=0.5, label='JULES')
plt.hist(d['classic'].data, color='orange', alpha=0.5, label='CLASSIC')
plt.hist(d['visit'].data, color='darkred', alpha=0.5, label='VISIT')
plt.hist(d['ssib4'].data, color='yellowgreen', alpha=0.5, label='SSIB4')
plt.xlabel('Burnt Area (Mha)')
plt.ylabel('Probability')
plt.legend()
plt.show()
'''




