#conda activate geopackage in terminal 
#python3
import geopandas as gpd
import matplotlib.pyplot as plt
import geoplot
from matplotlib.colors import TwoSlopeNorm
import pandas as pd


#1) Download the geopackage hexagons file from here: https://gis.stackexchange.com/questions/415139/world-hexmap-based-on-ipcc-ar6-matching-country-codes-with-wgi-reference-region
#hexagons = gpd.read_file('/scratch/cburton/scratch/ISIMIP3a/Data/zones.gpkg')

#2) Delete PAC and GIC in the a row
#Go to QGIS and open the file. Layer - Open Attribute Table. Select row and Delete. Can also add a new column here
#Save as '/scratch/cburton/scratch/ISIMIP3a/Data/WeightedZones.gpkg'


#3)Make a CSV file with the weights in Excel, as 2 columns (region abbreviation and Probability Ratio), with the regions in the same order as the geopackage file
#Save as "/scratch/cburton/scratch/ISIMIP3a/Data/RiskRatios_AR6regions.csv"


#4) To Make the Hexagon plot in the terminal:
hexagons = gpd.read_file('/scratch/cburton/scratch/ISIMIP3a/Data/WeightedZones.gpkg')
df = gpd.read_file("/scratch/cburton/scratch/ISIMIP3a/Data/RiskRatios_AR6regions.csv")
df['field_2'] = pd.to_numeric(df['field_2'])
hexagons['Model_results'] = df['field_2']# Add the Probability Ratio from df to the Model_results column in the hexagons gpkg


#Centre the colourbar
vmin, vmax, vcenter = hexagons['Model_results'].min(), hexagons['Model_results'].max(), 1.0
norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
cbar = plt.cm.ScalarMappable(norm=norm, cmap='RdBu_r',)
#Make the plot
ax = hexagons.plot()
hexagons.plot(column='Model_results', cmap='RdBu_r', norm=norm, linewidth=0.8, edgecolor='gray', legend=False, ax=ax)
hexagons.apply(lambda x: ax.annotate(text=x['label'], xy=x.geometry.centroid.coords[0], ha='center'), axis=1);
plt.colorbar(cbar,ticks=[0.8,1.0,1.5,2.0],orientation='horizontal', label='Probability Ratio')
plt.axis('off')
plt.show()






#Other things that might be useful

#Just the plot
hexagons = gpd.read_file('/scratch/cburton/scratch/ISIMIP3a/Data/WeightedZones.gpkg')
ax = hexagons.plot()
hexagons.plot(column='Model_results', ax=ax, cmap='RdBu_r', linewidth=0.8, edgecolor='gray', legend = False)
hexagons.apply(lambda x: ax.annotate(text=x['label'], xy=x.geometry.centroid.coords[0], ha='center'), axis=1);
plt.axis('off')
plt.show()

# Initial Set Up
df = gpd.read_file("/scratch/cburton/scratch/ISIMIP3a/Data/RiskRatios_AR6regions.csv")

hexagons['Model_results'] = df['field_2']
ax = hexagons.plot()
hexagons.plot(column='Model_results', ax=ax, cmap='YlOrRd', linewidth=0.8, edgecolor='gray', legend = True, legend_kwds={'loc': 'lower right'});
hexagons.apply(lambda x: ax.annotate(text=x['label'], xy=x.geometry.centroid.coords[0], ha='center'), axis=1);
plt.axis('off')
plt.show()


#Print all tables
cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
table_names = cursor.fetchall()
print(table_names)


# print the column names
cursor.execute("SELECT * FROM WeightedZones")
rows = cursor.fetchall()
for description in cursor.description:
    print(description[0], end="\t")
print()


#Save changes
hexagons.to_file('/scratch/cburton/scratch/ISIMIP3a/Data/WeightedZones.gpkg', driver='GPKG')


print(hexagons.head())
print(hexagons.label())

#plot a map
hexagons.plot(cmap='magma')
plt.show()


#Add region labels
ax = hexagons.plot()
hexagons.apply(lambda x: ax.annotate(text=x['label'], xy=x.geometry.centroid.coords[0], ha='center'), axis=1);
plt.show()

geoplot.choropleth(hexagons, hue=df['field_2'], cmap='magma')


#Use ID column as colour
hexagons.plot(column='id', ax=ax, cmap='YlOrRd', linewidth=0.8, edgecolor='gray', legend = True)



#Some links about geopandas: 
https://geobgu.xyz/py/geopandas1.html
https://gis.stackexchange.com/questions/342855/reading-geopackage-geometries-in-python
https://geopandas.org/en/stable/docs/user_guide/mapping.html
https://geopandas.org/en/stable/gallery/plotting_with_geoplot.html
https://gis.stackexchange.com/questions/336437/colorizing-polygons-based-on-color-values-in-dataframe-column
https://stackoverflow.com/questions/39816790/relocating-legend-from-geopandas-plot

