"""Assign agriculture GDP to land use layers in Jamaica
"""
import sys
import os
import subprocess

import rasterio
import rioxarray
import pandas as pd
import geopandas as gpd
# gpd._compat.USE_PYGEOS = True
# gpd.options.use_pygeos = True
from shapely.geometry import Point
import shapely
import numpy as np
from utils import *
from tqdm import tqdm
tqdm.pandas()

epsg_jamaica = 3448

def main(config):
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']


    crop_folders = ["spam2010v2r0_global_val_prod_agg.geotiff",
                    "spam2010v2r0_global_prod.geotiff",
                    "spam2010v2r0_global_yield.geotiff"]
    crop_strings = ["spam2010V2r0_global_V_agg_",
                        "spam2010V2r0_global_P_",
                        "spam2010V2r0_global_Y_"]
    crop_outputs = ["production","tonnage","yield"]
    for i, (crop_path,crop_field,crop_layer) in enumerate(list(zip(crop_folders,crop_strings,crop_outputs))):
        crop_data_path = os.path.join(incoming_data_path,
                                    "agriculture_data",
                                    crop_path,
                                    "JAM")
        all_crops = []
        all_fields = []
        for file in os.listdir(crop_data_path):
            if file.endswith(".tif"):
                crop_file = file.replace(".tif","")
                field_name = crop_file.replace(crop_field,"")
                all_fields.append(field_name)
                crop_raster_in = os.path.join(crop_data_path,
                                        f"{crop_file}.tif")
                crop_raster_out = os.path.join(crop_data_path,
                                        f"{crop_file}_reproject.tif")
                raster_rewrite(crop_raster_in,crop_raster_out)

                outCSVName = os.path.join(crop_data_path, f"{crop_file}.csv")
                subprocess.run(["gdal2xyz.py", '-csv', crop_raster_out, outCSVName])

                # Load points and convert to geodataframe with coordinates
                load_points = pd.read_csv(outCSVName, header=None, names=[
                                          'x', 'y', field_name], index_col=None)
                # load_points = load_points[load_points['crop'] > 0]
                os.remove(outCSVName)
                os.remove(crop_raster_out)

                if len(all_crops) > 0:
                    all_crops = pd.merge(all_crops,load_points,how="left",on=['x','y'])
                else:
                    all_crops = load_points.copy()
                
                del load_points
                print ("* Done with",file)

        all_crops["geometry"] = [Point(xy) for xy in zip(all_crops.x, all_crops.y)]
        crop_points = gpd.GeoDataFrame(all_crops, crs=f"EPSG:{epsg_jamaica}", geometry="geometry")
        crop_points["crop_id"] = crop_points.index.values.tolist()
        del all_crops
        print (crop_points) 

        crop_areas = create_voronoi_layer(crop_points,
                            "crop_id",epsg=epsg_jamaica)

        crop_areas = gpd.GeoDataFrame(pd.merge(crop_areas,crop_points[["crop_id"] + all_fields],
                                    how="left",on=["crop_id"]),
                                geometry="geometry",crs=f"EPSG:{epsg_jamaica}")

        crop_points.to_file(os.path.join(processed_data_path,
                                "agriculture_data",
                                "spam_agriculture_outputs.gpkg"),
                            layer=f"{crop_layer}_value",driver="GPKG")
        crop_areas.to_file(os.path.join(processed_data_path,
                                "agriculture_data",
                                "spam_agriculture_outputs.gpkg"),
                            layer=f"{crop_layer}_areas",driver="GPKG")
        all_fields_values = []
        for f in all_fields:
            all_fields_values.append((f,crop_areas[crop_areas[f] > 0][f].sum()))
        pd.DataFrame(all_fields_values,
                    columns=[f"{crop_layer}_column","value"]).to_csv(os.path.join(processed_data_path,
                                "agriculture_data",
                                f"{crop_layer}_column_keys.csv"),index=False)
        del crop_areas, crop_points

    crop_yields = gpd.read_file(os.path.join(processed_data_path,
                                "agriculture_data",
                                "spam_agriculture_outputs.gpkg"),
                            layer=f"tonnage_areas")
    crop_yields = crop_yields.to_crs(epsg=epsg_jamaica)
    print (crop_yields)
    crop_details = pd.read_csv(os.path.join(processed_data_path,
                                "agriculture_data",
                                "crop_details.csv")
                            )
    print (crop_details)
    tech_type = ["A","I","R"]
    poultry_crops = ["maiz","ocer","pmil","smil","soyb","sunf","whea"]
    
    crop_yields["crop_tons"] = 0

    all_crop_columns = []
    all_sector_columns = []
    for crop in crop_details.itertuples():
        crop_columns = [f"{crop.name.upper()}_{t}" for t in tech_type]
        sector_columns = [f"{crop.sector_code}_{crop.subsector_code}_{t}" for t in tech_type]
        poultry_columns = [f"A_12_{t}" for t in tech_type]
        all_sector_columns += sector_columns + poultry_columns
        for i,(cr,sc,pc) in enumerate(list(zip(crop_columns,sector_columns,poultry_columns))):
            if sc not in crop_yields.columns.values.tolist():
                crop_yields[sc] = 0
            if pc not in crop_yields.columns.values.tolist():
                crop_yields[pc] = 0
            crop_yields[cr] = crop_yields.apply(lambda x:x[cr] if x[cr] > 0 else 0,axis=1)
            crop_yields[f"{cr}_prod"] = crop.value_usd_per_ton*crop_yields[cr]
            crop_yields["crop_tons"] += crop_yields[cr]
            crop_yields[sc] += crop_yields[f"{cr}_prod"]
            if crop.name in poultry_crops:
                crop_yields[pc] += crop_yields[f"{cr}_prod"]
            all_crop_columns.append(f"{cr}_prod")

        print ("* Done with",crop.name,crop.value_usd_per_ton)

    crop_yields["crop_tons_persqm"] = crop_yields["crop_tons"]/crop_yields["areas"]
    crop_yields = crop_yields[["crop_id","areas","crop_tons_persqm","geometry"] + list(set(all_sector_columns)) + all_crop_columns]
    write_results = False
    if write_results is True:
        crop_yields.to_file(os.path.join(processed_data_path,
                                    "agriculture_data",
                                    "spam_agriculture_outputs.gpkg"),
                                layer=f"values_per_sqm",driver="GPKG")

    financial_year = 2019
    economic_output = pd.read_excel(os.path.join(processed_data_path,
                                        "macroeconomic_data",
                                        "detailed_sector_GVA_GDP_current_prices.xlsx"),
                        sheet_name="2019")
    economic_output.columns = [str(c).strip() for c in economic_output.columns.values.tolist()]
    economic_output["subsector_code"] = economic_output["subsector_code"].apply(str)
    totat_gva = economic_output[economic_output["sector_code"] == "GVA"][f"{financial_year}"].sum()
    total_tax =  economic_output[economic_output["sector_code"] == "TAX"][f"{financial_year}"].sum()
    tax_rate = 1.0*total_tax/totat_gva

    economic_output = economic_output[economic_output["sector_code"] == "A"]

    agri_land_use = gpd.read_file(os.path.join(
                        processed_data_path,
                        "land_type_and_use",
                        "jamaica_land_use_combined_with_sectors.gpkg"
                            ),
                        layer="areas"
                    )
    agri_land_use = agri_land_use.to_crs(epsg=epsg_jamaica)
    agri_land_use["land_id"] = agri_land_use.index.values.tolist()
    agri_land_use["land_id"] = agri_land_use.progress_apply(lambda x:f"land_{x.land_id}",axis=1)
    non_agri_land_use = agri_land_use[~((agri_land_use["sector_code_forest"] == "A") | (agri_land_use["sector_code_tnc"] == "A"))]
    non_agri_land_use["A_GDP"] = 0
    agri_land_use = agri_land_use[(agri_land_use["sector_code_forest"] == "A") | (agri_land_use["sector_code_tnc"] == "A")]
    agri_land_use["known_forest"] = agri_land_use.progress_apply(
                        lambda x:1 if (str(x.subsector_code_forest) == "20" and str(x.subsector_code_tnc) == "20") else 0,
                        axis=1)
    agri_forest = agri_land_use[agri_land_use["known_forest"] == 1]
    agri_land_use = agri_land_use[agri_land_use["known_forest"] == 0]
    for del_col in ["index","index_left","index_right"]:
        if del_col in agri_land_use.columns.values.tolist():
            agri_land_use.drop(del_col,axis=1,inplace=True)

    
    agri_areas = gpd.sjoin(agri_land_use,crop_yields,how="inner",predicate='intersects').reset_index()
    agri_areas.rename(columns={"geometry":"agri_geometry"},inplace=True)
    agri_areas = pd.merge(agri_areas, crop_yields[['crop_id','geometry']],how="left",on=["crop_id"])
    agri_areas["geom"] = agri_areas.progress_apply(lambda x:x["agri_geometry"].intersection(x["geometry"].buffer(0)),axis=1)
    agri_areas.drop(["agri_geometry","geometry"],axis=1,inplace=True)
    agri_areas.rename(columns={"geom":"geometry"},inplace=True)
    values_columns = ['land_id','forest_id','Classify', 
                        'LU_CODE','tnc_id', 
                        'NAME', 'TNCCODE', 
                        'global_id', 'global_LU_type',
                        "sector_code_forest",
                        "sector_code_tnc",'crop_tons_persqm'] + list(set(all_sector_columns)) + all_crop_columns + ["geometry"]
    agri_areas = agri_areas[values_columns]
    agri_areas["area_sqm"] = agri_areas.progress_apply(lambda x:x.geometry.area,axis=1)
    agri_areas["A_GDP"] = 0
    tot_gpd = 0
    for i,econ in economic_output.iterrows():
        econ_subsector_codes = str(econ["subsector_code"]).split(",")
        econ_codes = [f"A_{e}_A" for e in econ_subsector_codes if f"A_{e}_A" in agri_areas.columns.values.tolist()]
        if len(econ_codes) > 0:
            output = (1+tax_rate)*(1.0e6/365.0)*econ[f"{financial_year}"]
            tot_gpd += output
            agri_areas["weight"] = agri_areas[econ_codes].sum(axis=1)*agri_areas["area_sqm"]
            agri_areas["A_GDP"] += output*agri_areas["weight"]/agri_areas["weight"].sum()
    
    post_harvest_output = (1+tax_rate)*(1.0e6/365.0)*economic_output[
                        (economic_output["sector_code"] == "A"
                            ) & (economic_output["subsector_code"] == "14"
                            )][f"{financial_year}"].sum()
    tot_gpd += post_harvest_output
    print ("Sector A-14 output", post_harvest_output)
    agri_areas["crop_tons"] = agri_areas["crop_tons_persqm"]*agri_areas["area_sqm"]
    agri_areas["A_GDP"] += post_harvest_output*agri_areas["crop_tons"]/agri_areas["crop_tons"].sum()

    agri_areas = agri_areas.groupby(["land_id"])["A_GDP","crop_tons"].sum().reset_index()
    agri_areas = pd.merge(agri_land_use,agri_areas,how="left",on=["land_id"])


    forest_output = (1+tax_rate)*(1.0e6/365.0)*economic_output[
                        (economic_output["sector_code"] == "A"
                            ) & (economic_output["subsector_code"] == "20"
                            )][f"{financial_year}"].sum()
    tot_gpd += forest_output
    agri_forest["A_GDP"] = forest_output*agri_forest["area_m2"]/agri_forest["area_m2"].sum()
    values_columns = ['land_id','forest_id','Classify', 
                        'LU_CODE','tnc_id', 
                        'NAME', 'TNCCODE', 
                        'global_id', 'global_LU_type',"sector_code_forest",
                        "sector_code_tnc",'area_m2','A_GDP','geometry']
    agri_areas = pd.concat([agri_areas[values_columns + ['crop_tons']],
                            agri_forest[values_columns],
                            non_agri_land_use[values_columns]],
                            axis=0,ignore_index=True)
    agri_areas["crop_tons"] = agri_areas["crop_tons"].fillna(0)
    agri_areas["A_GDP"] = agri_areas["A_GDP"].fillna(0)
    
    agri_areas['GDP_persqm'] = agri_areas["A_GDP"]/agri_areas["area_m2"]
    agri_areas['crop_tons_persqm'] = agri_areas["crop_tons"]/agri_areas["area_m2"]
    agri_areas["GDP_unit"] = "JD/day"

    print ("* Given GDP",tot_gpd)
    print ("* Estimated GDP",agri_areas["A_GDP"].sum())

    agri_areas = gpd.GeoDataFrame(agri_areas,geometry="geometry",crs=f"EPSG:{epsg_jamaica}")
    agri_areas = remove_geometry_collections(agri_areas)
    write_results = False
    if write_results is True:
        agri_areas.to_file(os.path.join(processed_data_path,
                                    "agriculture_data",
                                    "agricuture_gdp.gpkg"),
                                layer="areas",driver="GPKG")
    del crop_yields, agri_forest,agri_land_use,economic_output

    tot_gpd = agri_areas["A_GDP"].sum()
    tot_area = agri_areas["area_m2"].sum()
    print ("* Estimated GDP",tot_gpd)
    print ("* Estimated Areas",tot_area)
    sector_df = gpd.read_file(os.path.join(incoming_data_path,
                                    "buildings",
                                    'buildings_assigned_economic_sectors_intermediate.gpkg'),
                            layer="commercial_sectors")
    sector_df = sector_df.to_crs(epsg=epsg_jamaica)
    sector_columns = [c for c in sector_df.columns.values.tolist() if "A_" in c[:2]]
    sector_df["A"] = sector_df[sector_columns].sum(axis=1)
    sector_df["A"] = sector_df.progress_apply(lambda x:1 if x["A"] > 0 else 0,axis=1)
    agri_buildings = gpd.sjoin(agri_areas[["land_id","GDP_persqm","geometry"]],
                                sector_df[sector_df["A"] == 1][["osm_id","geometry"]],
                            how="inner",predicate='intersects').reset_index()
    
    agri_buildings.rename(columns={"geometry":"landuse_geometry"},inplace=True)
    agri_buildings = pd.merge(agri_buildings, sector_df[['osm_id','geometry']],how="left",on=["osm_id"])
    agri_buildings["area_sqm"] = agri_buildings.progress_apply(lambda x:x["landuse_geometry"].intersection(x["geometry"].buffer(0)).area,axis=1)
    print ("* Building areas",agri_buildings["area_sqm"].sum())
    print ("* Building areas ratio",agri_buildings["area_sqm"].sum()/tot_area)
    agri_buildings["A_GDP_building"] = agri_buildings["GDP_persqm"]*agri_buildings["area_sqm"]
    agri_buildings_areas = agri_buildings.groupby(["land_id"])["A_GDP_building"].sum().reset_index()
    agri_buildings = agri_buildings.groupby(["osm_id"])["A_GDP_building"].sum().reset_index()
    print ("* GDP to buildings",agri_buildings["A_GDP_building"].sum())
    print ("* GDP ratio to buildings",agri_buildings["A_GDP_building"].sum()/tot_gpd)

    agri_areas = pd.merge(agri_areas,agri_buildings_areas,how="left",on=["land_id"])
    agri_areas["A_GDP_building"] = agri_areas["A_GDP_building"].fillna(0)
    agri_areas["GDP_building_ratio"] = agri_areas.progress_apply(lambda x:x["A_GDP_building"]/x["A_GDP"] if x["A_GDP"] > 0 else 0, axis=1)

    gpd.GeoDataFrame(agri_areas,
                    geometry="geometry",
                    crs=f"EPSG:{epsg_jamaica}").to_file(os.path.join(processed_data_path,
                                "agriculture_data",
                                "agricuture_gdp.gpkg"),
                            layer="areas",river="GPKG")
    agri_buildings.to_csv(os.path.join(processed_data_path,
                                    "agriculture_data",
                                    "building_agricuture_gdp.csv"),
                            index=False)

if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)