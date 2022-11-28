import iris
import matplotlib as mpl
import numpy as np
import iris.analysis.cartography
import h5py
from iris.coord_systems import GeogCS
import cartopy.crs as ccrs
import cf_units
import cftime
import iris.coord_categorisation
import matplotlib
import os
import glob
import sys
import pylab
import iris.plot as iplt
import iris.quickplot as qplt
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

cube = iris.load_cube('/scratch/cburton/scratch/ISIMIP3a/Data/jules*obsclim*.nc')
cube = UpdateTime(cube)
cube = ConstrainTime(cube)
print (cube)
GFED = iris.load_cube('/data/cr1/cburton/GFED/GFED4.1s/GFED4.1s_1997-2016.nc', iris.Constraint(cube_func=lambda x: x.var_name == 'basis_regions')).collapsed('time', iris.analysis.MEAN)
qplt.pcolormesh(GFED)
plt.show()
GFED = GFED.regrid(cube, iris.analysis.Linear())

template = cube.copy()
for i in range(len(template.coord('time').points)):
    template.data[i, :, :] = GFED.data[:, :]

template.units = cf_units.Unit(1)
template.rename('Basis_Regions')    

iris.save(template, '/scratch/cburton/scratch/ISIMIP3a/Data/GFED_Regions.nc')
print ('saved')
exit()




