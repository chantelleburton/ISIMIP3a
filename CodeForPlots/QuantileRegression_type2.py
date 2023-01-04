import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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
import iris.quickplot as qplt
from scipy import stats
import numpy.ma as ma
import matplotlib.colors as mcols
import statsmodels.api as sm
from scipy.stats import norm, uniform
import cartopy.feature as cfeature
from collections import OrderedDict 
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf


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



GFED = iris.load_cube('/scratch/cburton/scratch/ISIMIP3a/Data/GFED_Regions.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'basis_regions'))
folder = '/scratch/cburton/scratch/ISIMIP3a/Data/'
AllData = ['jules', 'classic','visit', 'ssib4', 'LPJ-GUESS-SPITFIRE','GFED4.1s', 'FireCCI5.1', 'FireCCILT11','counterclim/*jules', 'counterclim/*classic','counterclim/*visit', 'counterclim/*ssib4', 'counterclim/*LPJ-GUESS-SPITFIRE']
d = {}
time = np.arange(0,192) 
quantiles = np.arange(0.05, 0.96, 0.1)
df = pd.DataFrame({'q': quantiles})

for key in file_dict:
    for Data in AllData:
        cube = iris.load_cube(folder+Data+'*.nc')
        if Data == 'LPJ-GUESS-SPITFIRE' or Data == 'LPJ-GUESS-SIMFIRE' or Data == 'counterclim/*LPJ-GUESS-SPITFIRE':
            cube.data = ma.masked_where(np.isnan(cube.data),cube.data)
        elif Data == 'GFED4.1s' or Data == 'FireCCI5.1' or Data == 'FireCCILT11':
            cube = cube
        else:
          cube = UpdateTime(cube)
        cube = ConstrainTime(cube)
        if Data == 'classic' or Data == 'counterclim/*classic':
            cube = cube*100 #Convert frac to percent
        if Data == 'ssib4' or Data == 'counterclim/*ssib4':
            cube = cube*30 #Convert %/day to %/month
        cube.data = np.ma.masked_where(GFED.data != file_dict[key],cube.data)
        cube = CollapseToTimeseries(cube)
        cube = MakeAnomaly(cube).data

        BA = cube
        d = pd.DataFrame({'time': time, 'BA': BA})
        ols = smf.ols("BA ~ time", d).fit()
        ols = dict(a=ols.params["Intercept"], b=ols.params["time"])

        mod = smf.quantreg('BA ~ time', d)
        res = mod.fit(q=0.5)

        def fit_model(q):
            res = mod.fit(q=q)
            return res.params["Intercept"]

        models = [fit_model(x) for x in quantiles]
        df[Data] = models
        df[Data] = df[Data].astype(float)
        print('done',Data)

    print(df)

    plt.plot(df['q'], df['jules'], color="lightblue", label='JULES')
    plt.plot(df['q'], df['classic'], color="coral", label="classic")
    plt.plot(df['q'], df['visit'], color="sandybrown", label="visit")
    plt.plot(df['q'], df['ssib4'], color="lightgreen", label="ssib4")
    plt.plot(df['q'], df['LPJ-GUESS-SPITFIRE'], color="plum", label="LPJ-G-SPITFIRE")
    plt.plot(df['q'], df['GFED4.1s'], color="black", label='GFED4.1s')
    plt.plot(df['q'], df['FireCCI5.1'], color="blue", label="FireCCI5.1")
    plt.plot(df['q'], df['FireCCILT11'], color="darkblue", label="FireCCILT11")

    plt.plot(df['q'], df['counterclim/*jules'], linestyle='--', color="lightblue", label='counterclim')
    plt.plot(df['q'], df['counterclim/*classic'], linestyle='--', color="coral")
    plt.plot(df['q'], df['counterclim/*visit'], linestyle='--', color="sandybrown")
    plt.plot(df['q'], df['counterclim/*ssib4'], linestyle='--', color="lightgreen")
    plt.plot(df['q'], df['counterclim/*LPJ-GUESS-SPITFIRE'], linestyle='--', color="plum")

    plt.plot(df['q'], [ols["a"]] * 10, color="red", label="Ordinary Least Squares")

    plt.ylabel("Burnt Area Trend")
    plt.xlabel("Quantiles of the distribution")
    plt.title(key)
    plt.legend(ncol=2)
    plt.savefig('/scratch/cburton/scratch/ISIMIP3a/December/QR_'+key+'_OLS.png')
    plt.close()
    print ('saved', key)

exit()





         

