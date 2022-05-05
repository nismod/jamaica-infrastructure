import os
import sys
import geopandas as gpd
import pandas as pd
from summarise_utils import *

def jamaica_currency_conversion():
    """Conversion from J$ to US$
    """
    return 0.0067
    
def jamaica_sector_attributes():    
    sector_attributes = [
                            {
                                "sector":"transport",
                                "sector_gpkg":"rail.gpkg",
                                "sector_label":"Railways",
                                "edge_layer":"edges",
                                "node_layer":"nodes",
                                "area_layer":None,
                                "edge_id_column":"edge_id",
                                "node_id_column":"node_id",
                                "area_id_column":None,
                                "edge_classify_column":"status",
                                "node_classify_column":"asset_type",
                                "area_classify_column":None,
                                "edge_damage_filter_column":"status",
                                "node_damage_filter_column":"status",
                                "area_damage_filter_column":None,
                                "edge_damage_categories":["Functional"],
                                "node_damage_categories":["Functional"],
                                "area_damage_categories":None,
                                "edge_categories":["Functional","Non-Functional"],
                                "node_categories":["station"],
                                "area_categories":None,
                            },
                            {
                                "sector":"transport",
                                "sector_gpkg":"port_polygon.gpkg",
                                "sector_label":"Ports",
                                "edge_layer":None,
                                "node_layer":None,
                                "edge_id_column":None,
                                "node_id_column":None,
                                "area_id_column":"node_id",
                                "area_layer":"areas",
                                "edge_damage_filter_column":None,
                                "node_damage_filter_column":None,
                                "area_damage_filter_column":None,
                                "edge_damage_categories":None,
                                "node_damage_categories":None,
                                "area_damage_categories":None,
                                "edge_classify_column":None,
                                "node_classify_column":None,
                                "area_classify_column":"category",
                                "edge_categories":None,
                                "node_categories":None,
                                "area_categories":["transport"],
                            },
                            {
                                "sector":"transport",
                                "sector_gpkg":"airport_polygon.gpkg",
                                "sector_label":"Airports",
                                "edge_layer":None,
                                "node_layer":None,
                                "area_layer":"areas",
                                "edge_id_column":None,
                                "node_id_column":None,
                                "area_id_column":"node_id",
                                "edge_damage_filter_column":None,
                                "node_damage_filter_column":None,
                                "area_damage_filter_column":None,
                                "edge_damage_categories":None,
                                "node_damage_categories":None,
                                "area_damage_categories":None,
                                "edge_classify_column":None,
                                "node_classify_column":None,
                                "area_classify_column":"category",
                                "edge_categories":None,
                                "node_categories":None,
                                "area_categories":["transport"],
                            },
                            {
                                "sector":"transport",
                                "sector_gpkg":"roads.gpkg",
                                "sector_label":"Roads",
                                "edge_layer":"edges",
                                "node_layer":"nodes",
                                "area_layer":None,
                                "edge_id_column":"edge_id",
                                "node_id_column":"id",
                                "area_id_column":None,
                                "edge_classify_column":"road_class",
                                "node_classify_column":"asset_type",
                                "area_classify_column":None,
                                "edge_damage_filter_column":None,
                                "node_damage_filter_column":None,
                                "area_damage_filter_column":None,
                                "edge_damage_categories":None,
                                "node_damage_categories":None,
                                "area_damage_categories":None,
                                "edge_categories":["CLASS A","CLASS B","CLASS C","METRO","TRACK","OTHER"],
                                "node_categories":["bridge"],
                                "area_categories":None,
                            },
                            {
                                "sector":"energy",
                                "sector_gpkg":"electricity_network_v3.1.gpkg",
                                "sector_label":"Energy",
                                "edge_layer":"edges",
                                "node_layer":"nodes",
                                "area_layer":None,
                                "edge_id_column":"id",
                                "node_id_column":"id",
                                "area_id_column":None,
                                "edge_classify_column":"voltage_kV",
                                "node_classify_column":"subtype",
                                "area_classify_column":None,
                                "edge_damage_filter_column":None,
                                "node_damage_filter_column":None,
                                "area_damage_filter_column":None,
                                "edge_damage_categories":None,
                                "node_damage_categories":None,
                                "area_damage_categories":None,
                                "edge_categories":[138,69,24,12],
                                "node_categories":["diesel","gas","hydro",
                                                "substation","pole","demand"],
                                "area_categories":None,
                            },
                            {
                                "sector":"water",
                                "sector_gpkg":"", # We will integrate the NWC layers together for this
                                "sector_label":"Potable water",
                                "edge_layer":"edges",
                                "node_layer":"nodes",
                                "area_layer":None,
                                "edge_id_column":"edge_id",
                                "node_id_column":"node_id",
                                "area_id_column":None,
                                "edge_classify_column":"asset",
                                "node_classify_column":"asset_type",
                                "area_classify_column":None,
                                "edge_damage_filter_column":None,
                                "node_damage_filter_column":None,
                                "area_damage_filter_column":None,
                                "edge_damage_categories":None,
                                "node_damage_categories":None,
                                "area_damage_categories":None,
                                "edge_categories":["pipeline"],
                                "node_categories":["Booster Station","Catchment","Entombment",
                                                    "Filter Plant","Intake",
                                                    "Production Well","Pump Station",
                                                    "Relift Station","Reservoir",
                                                    "River Source","Spring","Storage Tank",
                                                    "Sump","Treatment Plant"],
                                "area_categories":None,
                            },
                            {
                                "sector":"water",
                                "sector_gpkg":"irrigation_assets_NIC.gpkg",
                                "sector_label":"Irrigation",
                                "edge_layer":"edges",
                                "node_layer":"nodes",
                                "area_layer":None,
                                "edge_id_column":"edge_id",
                                "node_id_column":"node_id",
                                "area_id_column":None,
                                "edge_classify_column":"asset_type",
                                "node_classify_column":"asset_type",
                                "area_classify_column":None,
                                "edge_damage_filter_column":None,
                                "node_damage_filter_column":None,
                                "area_damage_filter_column":None,
                                "edge_damage_categories":None,
                                "node_damage_categories":None,
                                "area_damage_categories":None,
                                "edge_categories":["canal","pipeline"],
                                "node_categories":["well"],
                                "area_categories":None,
                            },
                            {
                                "sector":"water",
                                "sector_gpkg":"waste_water_facilities_NWC.gpkg",
                                "sector_label":"Wastewater Treatment",
                                "edge_layer":None,
                                "node_layer":"nodes",
                                "area_layer":None,
                                "edge_id_column":None,
                                "node_id_column":"node_id",
                                "area_id_column":None,
                                "edge_classify_column":None,
                                "node_classify_column":"asset_type",
                                "area_classify_column":None,
                                "edge_damage_filter_column":None,
                                "node_damage_filter_column":None,
                                "area_damage_filter_column":None,
                                "edge_damage_categories":None,
                                "node_damage_categories":None,
                                "area_damage_categories":None,
                                "edge_categories":None,
                                "node_categories":["Sump","WW Pump Station","WW Relift Station","WW Treatment Plant"],
                                "area_categories":None,
                            },

                        ]

    return sector_attributes

def jamaica_port_and_airport_nodes():
    config = load_config()
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    figures_data_path = config['paths']['figures']


    transport_data_path = os.path.join(processed_data_path,
                        "networks",
                        "transport")
    ports = gpd.read_file(os.path.join(transport_data_path,"port_polygon.gpkg"),layer="areas").to_crs(epsg=JAMAICA_GRID_EPSG)
    ports["asset_type"] = "port"
    airports = gpd.read_file(os.path.join(transport_data_path,"airport_polygon.gpkg"),layer="areas").to_crs(epsg=JAMAICA_GRID_EPSG)
    airports["asset_type"] = "airport"
    nodes = pd.concat([ports[["node_id","name","asset_type","geometry"]],
                        airports[["node_id","name","asset_type","geometry"]]],
                        axis=0,ignore_index=True)
    
    nodes["centroid"] = nodes.geometry.centroid
    nodes.drop("geometry",axis=1,inplace=True)
    nodes.rename(columns={"centroid":"geometry"},inplace=True)
    nodes = gpd.GeoDataFrame(nodes,geometry="geometry",crs=f"EPSG:{JAMAICA_GRID_EPSG}")

    return nodes
