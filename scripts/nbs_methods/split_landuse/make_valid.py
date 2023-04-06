import geopandas as gpd
df = gpd.read_file("2013_landuse_Landcover.shp")
geom = df.geometry.buffer(0)
df.geometry = geom
df.to_file("2013_landuse_landcover.gpkg", driver="GPKG")

