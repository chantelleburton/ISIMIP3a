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
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf



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

GFED = iris.load_cube('/scratch/cburton/scratch/ISIMIP3a/Data/GFED_Regions.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'basis_regions'))
folder = '/scratch/cburton/scratch/ISIMIP3a/Data/'
models = ['jules']
d = {}
#for key in file_dict:
for model in models:
        cube = iris.load_cube(folder+model+'*_gswp3-w5e5_obsclim_histsoc_default_burntarea-total_global_monthly_1901_2019.nc')
        if model == 'LPJ-GUESS-SPITFIRE' or model == 'LPJ-GUESS-SIMFIRE':
            cube.data = ma.masked_where(np.isnan(cube.data),cube.data)
        else:
          cube = UpdateTime(cube)
        cube = ConstrainTime(cube)
        if model == 'classic':
            cube = cube*100 #Convert frac to percent
        if model == 'ssib4':
            cube = cube*30 #Convert %/day to %/month
        #cube.data = np.ma.masked_where(GFED.data != file_dict[key],cube.data)
        cube = CollapseToTimeseries(cube)
        #cube = MakeAnomaly(cube)


time = np.arange(0,192) 
BA = cube.data
df = pd.DataFrame({'time': time, 'BA': BA})
print(df.head())


mod = smf.quantreg('BA ~ time', df)
res = mod.fit(q=0.5)
print(res.summary())

quantiles = np.arange(0.05, 0.96, 0.1)
def fit_model(q):
    res = mod.fit(q=q)
    return [q, res.params["Intercept"], res.params["time"]] + res.conf_int().loc[
        "time"
    ].tolist()

models = [fit_model(x) for x in quantiles]
models = pd.DataFrame(models, columns=["q", "a", "b", "lb", "ub"])
ols = smf.ols("BA ~ time", df).fit()
ols_ci = ols.conf_int().loc["time"].tolist()
ols = dict(
    a=ols.params["Intercept"], b=ols.params["time"], lb=ols_ci[0], ub=ols_ci[1]
)
print(models)
print(ols)


### Plot type (1) - quantiles as horizontal lines
x = np.arange(time.min(), time.max())
get_y = lambda a, b: a + b * x
fig, ax = plt.subplots(figsize=(8, 6))
for i in range(models.shape[0]):
    y = get_y(models.a[i], models.b[i])
    ax.plot(x, y, linestyle="dotted", color="grey")

ax.plot(x, y, linestyle="dotted", color="grey", label='Quantiles (0.05-0.95)')
y = get_y(ols["a"], ols["b"])

ax.plot(x, y, color="red", label="Ordinary Least Squares")
ax.scatter(time, BA, alpha=0.2, label='Model data')
legend = ax.legend()
ax.set_xlabel("Time (Months since 2001)", fontsize=16)
ax.set_ylabel("BA z-score", fontsize=16)
plt.title('Global')
plt.show()



n = models.shape[0]
print("n=....",n)
exit()



         

