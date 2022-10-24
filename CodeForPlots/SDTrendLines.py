import iris
import matplotlib as mpl
import numpy as np
import iris.plot as iplt
import matplotlib
import matplotlib.pyplot as plt
from ascend import shape
import ascend
import cartopy.feature as cfeature
import numpy.ma as ma

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

##Change region here
region = 'SSA'

dataseries = ['jules','classic','visit','ssib4', 'LPJ-GUESS-SPITFIRE','GFED4.1s','FireCCI5.1','FireCCILT11']
d = {}
for data in dataseries:
    trend = iris.load_cube('/scratch/cburton/scratch/ISIMIP3a/SEPTEMBER_UPDATE/'+data+'_TrendCube.nc')
    trend = MaskRegion(trend, region=region).data
    #trend = MaskOcean(trend).data
    d[(data)] = [ trend.mean()-2*trend.std(), trend.mean()+2*trend.std() ]

x=np.arange(7)
plt.fill_between(x,d['GFED4.1s'][0],d['GFED4.1s'][1], color='lightblue', label='GFED')
plt.fill_between(x,d['FireCCI5.1'][0],d['FireCCI5.1'][1],  color='grey', label='CCI')
plt.fill_between(x,d['FireCCILT11'][0],d['FireCCILT11'][1], color='grey')
plt.vlines(x=1, ymin=d['jules'][0], ymax=d['jules'][1], color='black')
plt.vlines(x=2, ymin=d['classic'][0], ymax=d['classic'][1], color='black')
plt.vlines(x=3, ymin=d['visit'][0], ymax=d['visit'][1], color='black')
plt.vlines(x=4, ymin=d['ssib4'][0], ymax=d['ssib4'][1], color='black')
plt.vlines(x=5, ymin=d['ssib4'][0], ymax=d['LPJ-GUESS-SPITFIRE'][1], color='black')
plt.ylabel('Trend')
models=('JULES', 'CLASSIC', 'VISIT', 'SSIB4', 'LPJG-SPITFIRE')
x_pos = np.arange(len(models))+1
plt.xticks(x_pos, models)
plt.legend()
plt.title('2001-2016 Trend, SSA')
plt.show()



         

