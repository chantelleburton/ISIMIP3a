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
from collections import OrderedDict 

matplotlib.use('Agg')
import matplotlib.pyplot as plt


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


#Start text document - only do this once, so comment out if running again
#f = open('/scratch/cburton/scratch/ISIMIP3a/October/Attribution_PDF.dat','a')
#np.savetxt(f,('Region','Mean','95Percentile'),newline=' ',fmt='  %s')
#f.close()

###Load GFED regions
file_dict = OrderedDict([
    ('BONA',1),
    ('TENA',2),
    ('CEAM',3),
    ('NHSA',4),
    ('SHSA',5),
    ('EURO',6),
    ('MIDE',7),
    ('NHAF',8),
    ('SHAF',9),
    ('BOAS',10),
    ('CEAS',11),
    ('SEAS',12),
    ('EQAS',13),
    ('AUST',14)])

GFED = iris.load_cube('/scratch/cburton/scratch/ISIMIP3a/Data/GFED_Regions.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'basis_regions'))

##### Load Model Data ##### 
folder = '/scratch/cburton/scratch/ISIMIP3a/Data/'
models = ['jules', 'classic','visit', 'ssib4', 'LPJ-GUESS-SPITFIRE']
exps = ['obsclim', 'counterclim']
d = {}
for key in file_dict:
    print(key)
    for exp in exps:
        for model in models:
            cube = iris.load_cube(folder+model+'*_gswp3-w5e5_'+exp+'*_burntarea-total_global_monthly_1901_2019.nc')
            if model == 'LPJ-GUESS-SPITFIRE' or model == 'LPJ-GUESS-SIMFIRE':
                cube.data = ma.masked_where(np.isnan(cube.data),cube.data)
            else:
                cube = UpdateTime(cube)
            cube = ConstrainTime(cube)
            if model == 'classic':
                cube = cube*100 #Convert frac to percent
            if model == 'ssib4':
                cube = cube*30 #Convert %/day to %/month
            cube.data = np.ma.masked_where(GFED.data != file_dict[key],cube.data)
            d[(model+exp+key)] = CollapseBySpace(cube)

    ### Make arrays
    obsclim = np.append([d['julesobsclim'+key].data],[d['classicobsclim'+key].data,d['visitobsclim'+key].data,d['ssib4obsclim'+key].data,d['LPJ-GUESS-SPITFIREobsclim'+key].data])
    counterclim = np.append([d['julescounterclim'+key].data],[d['classiccounterclim'+key].data,d['visitcounterclim'+key].data,d['ssib4counterclim'+key].data,d['LPJ-GUESS-SPITFIREcounterclim'+key].data])

    ## Calculate RR for MEAN
    counterclimMEAN = np.mean(counterclim)  
    obsclimMEAN = np.mean(obsclim)  
    ALL = (np.count_nonzero(obsclim.data > counterclimMEAN))
    NAT = (np.count_nonzero(counterclim.data > counterclimMEAN))
    RRmean = ALL/NAT

    ## Calculate RR for 95 percentile
    counterclim95 = np.percentile(counterclim, 95)  
    obsclim95 = np.percentile(obsclim, 95)  
    ALL = (np.count_nonzero(obsclim.data > counterclim95))
    NAT = (np.count_nonzero(counterclim.data > counterclim95))
    RR95 = ALL/NAT

    f = open('/scratch/cburton/scratch/ISIMIP3a/October/Attribution_PDF.dat','a')
    np.savetxt(f,(key,RRmean,RR95),newline=' ',fmt='  %s')
    f.write('\n')
    f.close()

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
    plt.savefig('/scratch/cburton/scratch/ISIMIP3a/October/Attribution_PDF_'+key+'.png')
    plt.close()











