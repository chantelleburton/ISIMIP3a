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
from collections import OrderedDict 

cube = iris.load_cube('/scratch/cburton/scratch/ISIMIP3a/SEPTEMBER_UPDATE/jules_TrendCube.nc')
var_constraint = iris.Constraint(cube_func=lambda x: x.var_name == 'basis_regions')
GFED = iris.load_cube('/data/cr1/cburton/GFED/GFED4.1s/GFED4.1s_1997-2016.nc', var_constraint).collapsed(['time'], iris.analysis.MEAN)
GFED = GFED.regrid(cube, iris.analysis.Linear()) 

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

for key in file_dict:
    print(key)
    dataseries = ['jules','classic','visit','ssib4', 'LPJ-GUESS-SPITFIRE','GFED4.1s','FireCCI5.1','FireCCILT11']
    d = {}
    for data in dataseries:
        trend = iris.load_cube('/scratch/cburton/scratch/ISIMIP3a/October/'+data+'_TrendCube_Monthly.nc')
        trend = np.ma.masked_where(GFED.data != file_dict[key],trend.data)
        d[(data)] = [ trend.mean()-2*trend.std(), trend.mean()+2*trend.std() ]
        print (d[(data)])

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
    plt.title('2001-2016 Trend (Monthly), '+key)
    plt.savefig('/scratch/cburton/scratch/ISIMIP3a/October/SDtrend_'+key+'_Monthly.png')
    plt.close()

print('done')

         

