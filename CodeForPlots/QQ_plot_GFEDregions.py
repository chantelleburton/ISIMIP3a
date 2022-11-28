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
import cartopy.feature as cfeature
from collections import OrderedDict 


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


#cube = iris.load_cube('/scratch/cburton/scratch/ISIMIP3a/Data/jules*obsclim*.nc')
#GFED = iris.load_cube('/data/cr1/cburton/GFED/GFED4.1s/GFED4.1s_1997-2016.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'basis_regions'))
#GFED = GFED.regrid(cube, iris.analysis.Linear())
#GFED = ConstrainTime(GFED)

GFED = iris.load_cube('/scratch/cburton/scratch/ISIMIP3a/Data/GFED_Regions.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'basis_regions'))

##### Load Model Data ##### 
folder = '/scratch/cburton/scratch/ISIMIP3a/Data/'
AllData = ['jules', 'classic','visit', 'ssib4', 'LPJ-GUESS-SPITFIRE','GFED4.1s', 'FireCCI5.1', 'FireCCILT11']
d = {}
for key in file_dict:
    print(key)
    for Data in AllData:
        cube = iris.load_cube(folder+Data+'*.nc')
        if Data == 'LPJ-GUESS-SPITFIRE' or Data == 'LPJ-GUESS-SIMFIRE':
            cube.data = ma.masked_where(np.isnan(cube.data),cube.data)
        elif Data == 'GFED4.1s' or Data == 'FireCCI5.1' or Data == 'FireCCILT11':
            cube = cube
        else:
            cube = UpdateTime(cube)
        #cube = MonthlyToAnnual(cube)
        cube = ConstrainTime(cube)
        if Data == 'classic':
            cube = cube*100 #Convert frac to percent
        if Data == 'ssib4':
            cube = cube*30 #Convert %/day to %/month
        cube.data = np.ma.masked_where(GFED.data != file_dict[key],cube.data)
        cube = CollapseToTimeseries(cube)
        d[(Data)] = MakeAnomaly(cube)

    # Quantile-quantile / Scatter plot
    plt.figure()
    plt.plot(np.sort(d['GFED4.1s'].data), np.sort(d['jules'].data), marker='x', color='lightblue', label='JULES v GFED')
    plt.plot(np.sort(d['GFED4.1s'].data), np.sort(d['classic'].data), marker='x', color='coral', label='CLASSIC v GFED')
    plt.plot(np.sort(d['GFED4.1s'].data), np.sort(d['visit'].data), marker='x', color='sandybrown', label='VISIT v GFED')
    plt.plot(np.sort(d['GFED4.1s'].data), np.sort(d['ssib4'].data), marker='x', color='lightgreen', label='SSIB4 v GFED')
    plt.plot(np.sort(d['GFED4.1s'].data), np.sort(d['LPJ-GUESS-SPITFIRE'].data), marker='x', color='plum', label='SPITFIRE v GFED')

    plt.plot(np.sort(d['FireCCI5.1'].data), np.sort(d['jules'].data), color='lightblue', label='CCI51')
    plt.plot(np.sort(d['FireCCI5.1'].data), np.sort(d['classic'].data), color='coral')
    plt.plot(np.sort(d['FireCCI5.1'].data), np.sort(d['visit'].data), color='sandybrown')
    plt.plot(np.sort(d['FireCCI5.1'].data), np.sort(d['ssib4'].data), color='lightgreen')
    plt.plot(np.sort(d['FireCCI5.1'].data), np.sort(d['LPJ-GUESS-SPITFIRE'].data), color='plum')

    plt.plot(np.sort(d['FireCCILT11'].data), np.sort(d['jules'].data), linestyle='dashed', color='lightblue', label='CCI11')
    plt.plot(np.sort(d['FireCCILT11'].data), np.sort(d['classic'].data), linestyle='dashed', color='coral')
    plt.plot(np.sort(d['FireCCILT11'].data), np.sort(d['visit'].data), linestyle='dashed', color='sandybrown')
    plt.plot(np.sort(d['FireCCILT11'].data), np.sort(d['ssib4'].data), linestyle='dashed', color='lightgreen')
    plt.plot(np.sort(d['FireCCILT11'].data), np.sort(d['LPJ-GUESS-SPITFIRE'].data), linestyle='dashed', color='plum')

    axis_min=-2.5
    axix_max=4.0
    plt.xlabel('Observed burned area Z-score')
    plt.ylabel('FireMIP burned area Z-score')
    plt.plot([axis_min,axix_max],[axis_min,axix_max],color="k")
    plt.ylim(axis_min,axix_max)
    plt.xlim(axis_min,axix_max)
    plt.legend()
    plt.title(key)
    plt.savefig('/scratch/cburton/scratch/ISIMIP3a/October/QQ_'+key+'_Monthly.png')
    plt.close()

print('done')




         

