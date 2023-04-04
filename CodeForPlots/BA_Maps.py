import iris
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import cf_units
import cftime
from iris.plot import pcolormesh
import iris.quickplot as qplt
import iris.plot as iplt
from iris.coord_systems import GeogCS
from iris.coords import DimCoord



folder = '/scratch/cburton/scratch/ISIMIP3a/Data/'


def UpdateTime(cube):
    timeco = cube.coord('time')
    assert timeco.units == cf_units.Unit('months since 1901-01-01', calendar='360_day')
    timeco.units = cf_units.Unit('days since 1901-01-01', calendar='360_day')
    timeco.points = timeco.points * 30.
    return cube

def ConstrainTime(cube):
    date = iris.Constraint(time=lambda cell: 2002 <= cell.point.year <= 2019) 
    cube = cube.extract(date)
    return cube

def MonthlyToAnnual(cube):
    iris.coord_categorisation.add_season_year(cube, 'time', name='year')
    cube = cube.aggregated_by(['year'],iris.analysis.SUM)
    cube = cube.collapsed(['time'], iris.analysis.MEAN)
    return cube


GFED = iris.load_cube('~/GitHub/ISIMIP3a/Observations/GFED500m_Burned_Percentage.nc')*12/100
CCI = iris.load_cube(folder+'FireCCI5.1_Burned_Percentage.nc')*12/100
CLASSIC = iris.load_cube(folder+'classic*obsclim*.nc')*12
VISIT = iris.load_cube(folder+'visit*obsclim*.nc')*12/100
SSIB4 = iris.load_cube(folder+'ssib4*obsclim*.nc')*30*12/100
JULES = iris.load_cube(folder+'jules*obsclim*.nc')*12/100
SPITFIRE = iris.load_cube(folder+'*SPITFIRE*obsclim*.nc')*12/100
SIMFIRE = iris.load_cube(folder+'*SIMFIRE*obsclim*.nc')*12/100


AllObs = (GFED, CCI)
for ob in AllObs:
    ob = ConstrainTime(ob)
    ob = MonthlyToAnnual(ob)

AllModels = (CLASSIC, VISIT, SSIB4, JULES, SPITFIRE, SIMFIRE)
for model in AllModels:
    model = UpdateTime(model)
    model = ConstrainTime(model)
    model = MonthlyToAnnual(model)


####### Make Plots ######
plt.subplot(3,3,1)
iplt.pcolormesh(GFED, vmin=0, vmax=0.4, cmap='Reds')
plt.gca().coastlines()
plt.title('GFED 500m')

plt.subplot(3,3,2)
iplt.pcolormesh(CCI, vmin=0, vmax=0.4, cmap='Reds')
plt.gca().coastlines()
plt.title('Fire CCI')

plt.subplot(3,3,4)
iplt.pcolormesh(VISIT, vmin=0, vmax=0.4, cmap='Reds')
plt.gca().coastlines()
plt.title('VISIT')

plt.subplot(3,3,5)
iplt.pcolormesh(SSIB4, vmin=0, vmax=0.4, cmap='Reds')
plt.gca().coastlines()
plt.title('SSIB4-TRIFFID')

plt.subplot(3,3,6)
iplt.pcolormesh(CLASSIC, vmin=0, vmax=0.4, cmap='Reds')
plt.gca().coastlines()
plt.title('CLASSIC')

plt.subplot(3,3,7)
iplt.pcolormesh(JULES, vmin=0, vmax=0.4, cmap='Reds')
plt.gca().coastlines()
plt.title('JULES')

plt.subplot(3,3,8)
iplt.pcolormesh(SPITFIRE,  vmin=0, vmax=0.4, cmap='Reds')
plt.gca().coastlines()
plt.title('LPJ-GUESS-SPITFIRE')

plt.subplot(3,3,9)
mesh = iplt.pcolormesh(SIMFIRE, vmin=0, vmax=0.4, cmap='Reds')
plt.gca().coastlines()
plt.title('LPJ-GUESS-SIMFIRE-BLAZE')

                         #Position: Left, Bottom, width, height
colorbar_axes = plt.gcf().add_axes([0.27, 0.05, 0.45, 0.04])
colorbar = plt.colorbar(mesh, colorbar_axes, orientation='horizontal', label='fraction')
plt.tight_layout()

plt.show()


