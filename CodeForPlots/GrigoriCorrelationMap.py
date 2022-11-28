import iris
import matplotlib as mpl
import numpy as np
from iris.plot import pcolormesh
import iris.analysis.cartography
import iris.plot as iplt
from iris.coord_systems import GeogCS
import cartopy.crs as ccrs
import cf_units
import iris.coord_categorisation
import matplotlib
import matplotlib.pyplot as plt
import iris.quickplot as qplt
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
from cftime import date2num
from cftime import num2date
import os.path
import numpy.ma as ma
import cartopy.io.shapereader as shpreader
from ascend import shape
import ascend
import matplotlib.colors as mcols
import cartopy.feature as cfeature


def MaskRegion(cube, region):
    regions = ascend.shape.load_shp(ascend.EXAMPLE_GIORGI, short_name=region)
    domain_shape = regions.unary_union()   
    MaskedCube = domain_shape.mask_cube(cube)
    return MaskedCube

def MaskOcean(cube):
    landfr_cube = iris.load_cube('/scratch/cburton/scratch/ISIMIP3a/Data/mask_latlon2d.nc', 'land_sea_mask')
    oceanmask = landfr_cube.copy()
    oceanmask.data = ma.masked_less(landfr_cube.data, 0.5)

    cs_new = iris.coord_systems.GeogCS(6371229.0)
    coords = ('longitude', 'latitude')
    dataseries = (cube, oceanmask)
    for coord in coords:
        for data in dataseries:
            data.coord(coord).coord_system = cs_new
            data.coord(coord).long_name=None

    Cube = cube*oceanmask
    return Cube                        

def MonthlyToAnnual(cube):
    iris.coord_categorisation.add_year(cube, 'time', name='year')
    cube = cube.aggregated_by(['year'],iris.analysis.SUM)
    return cube

region = 'SSA'


#Plotting
minimum_log_level = 0.02
maximum_scale_level = 2.0
cmap='RdBu_r'
anom_norm = mcols.SymLogNorm(
    linthresh=minimum_log_level,
    linscale=0.01,
    vmin=-maximum_scale_level,
    vmax=maximum_scale_level)

tick_levels = [-2, -0.2, 0.0, 0.2, 2]


##### Load Data ##### 

dataseries = ['LPJ-GUESS-SPITFIRE']
#dataseries = ['jules','classic','visit','ssib4','GFED','FireCCI51','FireCCILT11']
d = {}

for data in dataseries:
    trend = iris.load_cube('/scratch/cburton/scratch/ISIMIP3a/Data/TrendCubes/'+data+'_TrendCube.nc')
    #trend = MaskOcean(trend)
    d[(data)] = MaskRegion(trend, region=region)

    plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=-35))
    ax.set_extent([-30.0, -85.0, 15.0, -60.0])
    ax.coastlines(resolution='110m')
    iplt.pcolormesh(d[(data)], cmap=cmap, norm=anom_norm)
    bar = plt.colorbar(ticks=tick_levels, orientation="horizontal")
    bar.set_ticklabels(tick_levels)
    plt.title(data)
    plt.savefig('/scratch/cburton/scratch/ISIMIP3a/SEPTEMBER_UPDATE/'+data+'.png')

exit()








'''
dataseries = ['jules','classic','visit','ssib4','GFED','FireCCI51','FireCCILT11']
d = {}

for data in dataseries:
    trend = iris.load_cube('/scratch/cburton/scratch/ISIMIP3a/Data/'+data+'_TrendCube.nc')
    trend = MaskOcean(trend)
    d[(data)] = MaskRegion(trend, region=region)

    ax = plt.subplot(3,3,n)
    iplt.pcolormesh(d[data], cmap=cmap, norm=anom_norm)
    ax.set_extent([-90, 10, 5, 85], crs=ccrs.PlateCarree())
    ax.coastlines(lw=0.2)
    cbar = plt.colorbar(orientation='horizontal')
    plt.title(str(data))
    n = n+1

plt.show()  
exit()

plt.subplot(3,3,1)
ax.set_extent([-30.0, -85.0, 15.0, -60.0])
iplt.pcolormesh(d['GFED'], cmap=cmap, norm=anom_norm)
plt.gca().coastlines()
plt.colorbar(orientation='horizontal')
plt.title('GFED4s')

plt.subplot(3,3,2)
iplt.pcolormesh(d['FireCCI51'], cmap=cmap, norm=anom_norm)
plt.gca().coastlines()
plt.colorbar(orientation='horizontal')
plt.title('CCI51')

plt.subplot(3,3,3)
iplt.pcolormesh(d['FireCCILT11'], cmap=cmap, norm=anom_norm)
plt.gca().coastlines()
plt.colorbar(orientation='horizontal')
plt.title('CCILT11')

plt.subplot(3,3,4)
iplt.pcolormesh(d['jules'], cmap=cmap, norm=anom_norm)
plt.gca().coastlines()
plt.colorbar(orientation='horizontal')
plt.title('JULES')

plt.subplot(3,3,5)
iplt.pcolormesh(d['classic'], cmap=cmap, norm=anom_norm)
plt.gca().coastlines()
plt.colorbar(orientation='horizontal')
plt.title('CLASSIC')

plt.subplot(3,3,5)
iplt.pcolormesh(d['visit'], cmap=cmap, norm=anom_norm)
plt.gca().coastlines()
plt.colorbar(orientation='horizontal')
plt.title('VISIT')

plt.subplot(3,3,5)
iplt.pcolormesh(d['ssib4'], cmap=cmap, norm=anom_norm)
plt.gca().coastlines()
plt.colorbar(orientation='horizontal')
plt.title('ssib4')


plt.show()

'''







'''
plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=-35))
ax.set_extent([-30.0, -85.0, 15.0, -60.0])
ax.coastlines(resolution='110m')
qplt.pcolormesh(GFED)
ax.coastlines()
plt.show()
exit()
'''


