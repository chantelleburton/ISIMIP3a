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

#GFED = iris.load_cube('/scratch/cburton/scratch/ISIMIP3a/Data/GFED_Regions.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'basis_regions'))

cube = iris.load_cube('/scratch/cburton/scratch/ISIMIP3a/Data/jules*obsclim*.nc')
GFED = iris.load_cube('/data/cr1/cburton/GFED/GFED4.1s/GFED4.1s_1997-2016.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'basis_regions'))
GFED = GFED.regrid(cube, iris.analysis.Linear())
GFED = ConstrainTime(GFED)


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

##### Load Data ##### 
folder = '/scratch/cburton/scratch/ISIMIP3a/Data/'
AllData = ['jules','GFED4.1s', 'FireCCI5.1', 'FireCCILT11','classic','visit', 'ssib4', 'LPJ-GUESS-SPITFIRE']
colors= ('lightblue','black', 'blue', 'darkblue', 'coral', 'sandybrown', 'lightgreen', 'plum') 
d = {}
n=0
for key in file_dict:
    for Data in AllData:
        print(Data)
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
        cube.data = np.ma.masked_where(GFED.data != file_dict[key],cube.data)
        cube = CollapseBySpace(cube)
        d[(Data)] = MakeAnomaly(cube)
        sns.distplot(d[Data].data, hist=False, kde=True, 
                 bins=int(180/5), color = colors[n % len(colors)], 
                 hist_kws={'edgecolor':colors[n % len(colors)]},
                 kde_kws={'linewidth': 4, 'color': colors[n % len(colors)], 'label':Data})
        n=n+1
    plt.legend(ncol=1, loc='best')
    plt.xlabel('Normalised Burnt Area')
    plt.ylabel('Probability Density')
    plt.title(key)
    plt.savefig('/scratch/cburton/scratch/ISIMIP3a/October/PDF_'+key+'_Annual.png')
    plt.close()

print('done')

'''
colors= ('lightblue','black', 'blue', 'darkblue', 'coral', 'sandybrown', 'lightgreen', 'plum') 
n=0
for test in np.arange(8):
    color = colors[n % len(colors)]
    print (color)
    n=n+1

exit()
'''



