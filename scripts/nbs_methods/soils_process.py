"""Assign areas to different soil types for Jamaica
    Find the areas by different classes as totals and percentages
"""
import sys
import os
import pandas as pd
import geopandas as gpd
import numpy as np
from utils import *
from tqdm import tqdm
tqdm.pandas()

# The projection system for Jamaica
epsg_jamaica = 3448
m2_to_hectare = 0.0001

def get_areas(gdf):
    gdf = gdf.to_crs(epsg=epsg_jamaica)
    gdf["area_m2"] = gdf.progress_apply(lambda x:x.geometry.area,axis=1)
    gdf["area_hectares"] = m2_to_hectare*gdf["area_m2"]
    gdf["percentage"] = 100.0*gdf["area_m2"]/gdf["area_m2"].sum() 

    return gdf

def round_and_convert_to_thousands(gdf,columns_rounding):
    for i,(col_name,value_precision) in enumerate(columns_rounding):
        gdf[col_name] = gdf.apply(lambda x: "{:,}".format(round(x[col_name],value_precision)), axis=1)
    return gdf

def group_areas_and_find_totals(gdf,groupby_column,columns_rename,columns_rounding): 
    gdf = get_areas(gdf)
    gdf = gdf.groupby([groupby_column])["area_m2","area_hectares","percentage"].sum().reset_index()  

    gdf.columns = columns_rename
    gdf = round_and_convert_to_thousands(gdf,columns_rounding)
    return gdf

def main(config):
    # Specify the data folder paths where all datasets as stored
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']

    nbs_input_data_path = os.path.join(incoming_data_path,"nbs")
    
    soil_layer = gpd.read_file(os.path.join(nbs_input_data_path,"nsmdb-soils.gpkg"))
    soil_layer["classes"] = soil_layer.progress_apply(lambda x:str(x["TEXTURE1"]).lower(),axis=1)
    print (soil_layer)
    soil_classification = pd.read_excel(os.path.join(incoming_data_path,
                                                "nbs","NbS information.xlsx"),
                                sheet_name="Soils") 
    soil_classification["classes"] = soil_classification.progress_apply(lambda x:str(x["Type"]).lower(),axis=1)
    soil_layer = pd.merge(soil_layer,soil_classification,how="left",on=["classes"])
    soil_layer["Permeability class"] = soil_layer["Permeability class"].fillna("Unknown")
    soil_layer["Weighting"] = soil_layer["Weighting"].fillna(0)
    soil_layer.drop("classes",axis=1,inplace=True)

    soil_layer = gpd.GeoDataFrame(soil_layer,geometry="geometry",crs=f"EPSG:{epsg_jamaica}")
    soil_layer.to_file(os.path.join(processed_data_path,"nbs","nsmdb-soils.gpkg"),driver="GPKG")

    excel_writer = pd.ExcelWriter(os.path.join(processed_data_path,
                                "nbs",
                                "soil_classes_with_areas.xlsx"))

    soil_grouping_column = "TEXTURE1"
    soil_by_type = group_areas_and_find_totals(soil_layer,soil_grouping_column,
                                        ["Soil Type",
                                        "Area (m2)", 
                                        "Area (Hectares)",
                                        "Area (% of total mainland area)"],
                                        [("Area (m2)",0), 
                                        ("Area (Hectares)",0),
                                        ("Area (% of total mainland area)",2)])
    print (soil_by_type)

    soil_by_type.to_excel(excel_writer,"soil_type", index=False)

    soil_grouping_column = "Permeability class"
    soil_by_class = group_areas_and_find_totals(soil_layer,soil_grouping_column,
                                        ["Permeability class",
                                        "Area (m2)", 
                                        "Area (Hectares)",
                                        "Area (% of total mainland area)"],
                                        [("Area (m2)",0), 
                                        ("Area (Hectares)",0),
                                        ("Area (% of total mainland area)",2)])
    print (soil_by_class)

    soil_by_class.to_excel(excel_writer,"soil_class", index=False)
    excel_writer.save()



if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)