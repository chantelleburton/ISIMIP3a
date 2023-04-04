
# On the command line:
#conda create --name mamba --channel conda-forge mamba
#conda activate mamba
#mamba create --name myenvname --channel conda-forge regionmask cartopy pygeos
#conda deactivate
#conda activate myenvname


import cartopy.crs as ccrs
import geopandas as gp
import regionmask
import matplotlib.pyplot as plt



mask = regionmask.defined_regions.ar6.all
#https://regionmask.readthedocs.io/en/stable/defined_scientific.html
print(mask)

text_kws = dict(color="#67000d", fontsize=7, bbox=dict(pad=0.2, color="w"))
regionmask.defined_regions.ar6.all.plot(
    text_kws=text_kws, label_multipolygon="all")
plt.show()
exit()



         

