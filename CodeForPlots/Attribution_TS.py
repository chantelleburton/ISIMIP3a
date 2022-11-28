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

'''
x = np.arange(2001,2017)
test = [3,4,5,6,3,4,5,7,8,9, 8, 7, 4, 5, 6, 7]
test1 = [13,14,15,16,13,14,15,17,18,19, 18, 17, 14, 15, 16, 17]
test2 = [2,5,6,7,12,24,14,11,25,26,27,23,11,12,13,14]
test3 = [23,24,25,26,23,24,25,27,28,19, 18, 17, 14, 15, 16, 17]
test4 = [25,26,28,29,30,24,25,27,28,19, 18, 17, 14, 15, 16, 17]

AllData = [test, test1]
print(type(AllData))
AllData2 = [test2, test3]
         
def box_plot(data, fill_color):
    bp = plt.boxplot(data, whis=(0,100), widths=0, patch_artist=True)
    for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
        plt.setp(bp[element], color=fill_color)
    for patch in bp['boxes']:
        patch.set(facecolor=fill_color) 
    for cap in bp['caps']:
        cap.set(xdata=cap.get_xdata() + (-0.1,0.1))      
    return bp

NewObsClim = np.array([])
NewCounterClim = []
for n in np.arange(1,16):
    CounterClim  = [test[n],test1[n],test2[n],test3[n],test4[n]]
    NewCounterClim.append(CounterClim)

print(type(NewCounterClim))
box_plot(NewCounterClim, 'blue')
plt.show()
exit()
#box_plot(AllData, 'blue')#, whis=(0,100), widths=0)
#box_plot(AllData2, 'orange')#, whis=(0,100), widths=0)
y_pos = np.arange(1,3)
xs = np.arange(2001, 2003)
plt.xticks(y_pos, xs)
plt.show()
exit()
'''


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

folder = '/scratch/cburton/scratch/ISIMIP3a/Data/'
d = {}
x = np.arange(2001,2017)

##### Load Model Data ##### 
models = ['jules', 'classic','visit', 'ssib4', 'LPJ-GUESS-SPITFIRE', 'LPJ-GUESS-SIMFIRE']
Exps = ['obsclim','counterclim']
color='lightblue'
for Exp in Exps:
    for model in models:
        cube = iris.load_cube(folder+model+'*_gswp3-w5e5_'+Exp+'_histsoc_*_burntarea-total_global_monthly_1901_2019.nc')
        if model == 'LPJ-GUESS-SPITFIRE' or model == 'LPJ-GUESS-SIMFIRE':
            cube.data = ma.masked_where(np.isnan(cube.data),cube.data)
        else:
           cube = UpdateTime(cube)
        cube = MonthlyToAnnual(cube)
        cube = ConstrainTime(cube)
        if model == 'classic':
            cube = cube*100 #Convert frac to percent
        if model == 'ssib4':
            cube = cube*30 #Convert %/day to %/month
        cube = CollapseToTimeseries(cube) 
        d[(model+Exp)] = MakeAnomaly(cube)
        #plt.scatter(x,d[(model+Exp)].data, color=color)
    color='gold'


Obsclim = (d['classicobsclim'].data,d['julesobsclim'].data,d['visitobsclim'].data,d['ssib4obsclim'].data,d['LPJ-GUESS-SPITFIREobsclim'].data)
ObsclimMax=np.amax(Obsclim, axis=0)
ObsclimMin=np.amin(Obsclim, axis=0)
ObsclimMean=np.mean(Obsclim, axis=0)

Counterclim = (d['classiccounterclim'].data,d['julescounterclim'].data,d['visitcounterclim'].data,d['ssib4counterclim'].data,d['LPJ-GUESS-SPITFIREcounterclim'].data)
CFMax=np.amax(Counterclim, axis=0)
CFMin=np.amin(Counterclim, axis=0)
CFMean=np.mean(Counterclim, axis=0)

'''
##### Make scatter Plot ######
plt.scatter(x,ObsclimMean, color='blue', label='Hist')
plt.scatter(x,CFMean, color='darkorange',label='Counterclim')
plt.legend()
plt.ylabel('Z-score')
plt.title('Annual total burned area')
plt.show()
'''
'''
##### Make Range Plot ######
x = np.arange(2001,2017)
plt.fill_between(x, ObsclimMin, ObsclimMax, color='lightblue',alpha=0.5, label='Hist range')
plt.fill_between(x, CFMin, CFMax, color='orange', alpha=0.5, label='Counterclim range')
plt.legend()
plt.ylabel('Z-score')
plt.title('Annual total burned area Z-score')
plt.show()
'''

##### Make whikers Plot ######
def box_plot(data, fill_color):
    bp = plt.boxplot(data, whis=(0,100), widths=0, patch_artist=True)
    for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
        plt.setp(bp[element], color=fill_color)
    for patch in bp['boxes']:
        patch.set(facecolor=fill_color) 
    for cap in bp['caps']:
        cap.set(xdata=cap.get_xdata() + (-0.1,0.1))      
    return bp

NewObsClim = []
NewCounterClim = []

for n in np.arange(0,16):
    NewObsClim1 = [ d['classicobsclim'].data[n],d['julesobsclim'].data[n],d['visitobsclim'].data[n],d['ssib4obsclim'].data[n],d['LPJ-GUESS-SPITFIREobsclim'].data[n] ]
    NewObsClim.append(NewObsClim1)

    NewCounterClim1 = [ d['classiccounterclim'].data[n],d['julescounterclim'].data[n],d['visitcounterclim'].data[n],d['ssib4counterclim'].data[n],d['LPJ-GUESS-SPITFIREcounterclim'].data[n] ]
    NewCounterClim.append(NewCounterClim1)

box_plot(NewObsClim, 'blue')
box_plot(NewCounterClim, 'orange')
xf = np.arange(1,17)
plt.scatter(xf,ObsclimMean, color='blue', label='HIST')
plt.scatter(xf,CFMean, color='orange', label='Counterfact')
plt.xticks(xf,x)
plt.legend()
plt.ylabel('Z-score')
plt.title('Annual total burned area')
plt.show()

    


