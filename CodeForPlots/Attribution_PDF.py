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

matplotlib.use('Agg')
import matplotlib.pyplot as plt



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
    #cube = cube.collapsed(coords, iris.analysis.PERCENTILE, percent=[95])/1E12*100/100 # Convert to Mkm2 and then to Mha  
    return cube

def CollapseByTime(cube):
    cube = cube.collapsed('time', iris.analysis.MEAN) 
    #cube = cube.collapsed('time', iris.analysis.PERCENTILE, percent=[95]) 
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
models = ['LPJ-GUESS-SPITFIRE', 'jules', 'classic','visit', 'ssib4']
exps = ['obsclim', 'counterclim']
d = {}
for exp in exps:
    for model in models:
        cube = iris.load_cube(folder+model+'*_gswp3-w5e5_'+exp+'*_burntarea-total_global_monthly_1901_2019.nc')
        if model == 'LPJ-GUESS-SPITFIRE' or model == 'LPJ-GUESS-SIMFIRE':
            cube.data = ma.masked_where(np.isnan(cube.data),cube.data)
        else:
            cube = UpdateTime(cube)
        if model == 'classic':
           cube = cube*100 #Convert frac to percent
        elif model == 'ssib4':
            cube = cube*30 #Convert %/day to %/month
        cube = ConstrainTime(cube)
        #cube = MaskRegion(cube, region=region)
        d[(model+exp)] = CollapseBySpace(cube)

### Make arrays
obsclim = np.append([d['julesobsclim'].data],[d['classicobsclim'].data,d['visitobsclim'].data,d['ssib4obsclim'].data,d['LPJ-GUESS-SPITFIREobsclim'].data])
counterclim = np.append([d['julescounterclim'].data],[d['classiccounterclim'].data,d['visitcounterclim'].data,d['ssib4counterclim'].data,d['LPJ-GUESS-SPITFIREcounterclim'].data])

print (obsclim.shape)
print (len(counterclim.data))
print (len(obsclim.data))
print (len(d['julesobsclim'].data))


## Calculate RR for 95 percentile
counterclim95 = np.percentile(counterclim, 95)  
obsclim95 = np.percentile(obsclim, 95)  
print(counterclim95)

ALL = (np.count_nonzero(obsclim.data > counterclim95))
print (ALL)
NAT = (np.count_nonzero(counterclim.data > counterclim95))
print (NAT)
RR = ALL/NAT
print('95th percentile',RR)

## Calculate RR for MEAN
counterclimMEAN = np.mean(counterclim)  
obsclimMEAN = np.mean(obsclim)  
print(counterclimMEAN)

ALL = (np.count_nonzero(obsclim.data > counterclimMEAN))
print (ALL)
NAT = (np.count_nonzero(counterclim.data > counterclimMEAN))
print (NAT)

RR = ALL/NAT
print('mean',RR)



##Plot style 1
sns.distplot(obsclim, kde=False,fit=stats.genextreme, color="red",label='HIST',fit_kws={"linewidth":2.5,"color":"darkred"})
sns.distplot(counterclim, kde=False,fit=stats.genextreme, color="green",label='CF',fit_kws={"linewidth":2.5,"color":"darkgreen"})
plt.axvline(x=obsclim95, color='red',linestyle='--')
plt.axvline(x=counterclim95, color='green',linestyle='--')
plt.axvline(x=obsclimMEAN, color='red')
plt.axvline(x=counterclimMEAN, color='green')
plt.legend()
plt.xlabel('Burnt Area')
plt.ylabel('Probability Density')
plt.savefig('/home/h01/cburton/GitHub/ISIMIP3a/Plots/Attribution_PDF.png')
plt.close()
plt.show()










