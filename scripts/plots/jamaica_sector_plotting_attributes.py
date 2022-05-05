import os
import sys
import geopandas as gpd
import pandas as pd
from plot_utils import *

def jamaica_currency_conversion():
    """Conversion from J$ to US$
    """
    return 0.0068
    
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
                                "edge_categories_colors":["#238b45","#969696"],
                                "node_categories_colors":["#00441b"],
                                "area_categories_colors":None,
                                "edge_categories_labels":["Functional","Non-Functional"],
                                "node_categories_labels":["Station"],
                                "area_categories_labels":None,
                                "edge_categories_linewidth":[1.0,1.0],
                                "edge_categories_zorder":[10,10],
                                "node_categories_markersize":[10.0],
                                "node_categories_marker":["."],
                                "node_categories_zorder":[11],
                                "node_loss_result": "rail_stations",
                                "edge_loss_result":"road_rail_edges"
                            },
                            {
                                "sector":"transport",
                                "sector_gpkg":"port_polygon.gpkg",
                                "sector_label":"Ports",
                                "edge_layer":None,
                                "node_layer":None,
                                "edge_id_column":None,
                                "node_id_column":"node_id",
                                "area_id_column":"node_id",
                                "area_layer":"areas",
                                "edge_damage_filter_column":None,
                                "node_damage_filter_column":None,
                                "area_damage_filter_column":None,
                                "edge_damage_categories":None,
                                "node_damage_categories":None,
                                "area_damage_categories":None,
                                "edge_classify_column":None,
                                "node_classify_column":"category",
                                "area_classify_column":"category",
                                "edge_categories":None,
                                "node_categories":None,
                                "area_categories":["transport"],
                                "edge_categories_colors":None,
                                "node_categories_colors":None,
                                "area_categories_colors":["#08306b"],
                                "edge_categories_labels":None,
                                "node_categories_labels":None,
                                "area_categories_labels":["PORT AREAS"],
                                "edge_categories_linewidth":None,
                                "edge_categories_zorder":None,
                                "node_categories_markersize":[15.0],
                                "node_categories_marker":["o"],
                                "node_categories_zorder":[11],
                                "node_loss_result": "ports",
                                "edge_loss_result":None,
                                "node_loss_result":"ports",
                                "edge_loss_result":None
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
                                "edge_categories_colors":None,
                                "node_categories_colors":None,
                                "area_categories_colors":["#8c510a"],
                                "edge_categories_labels":None,
                                "node_categories_labels":None,
                                "area_categories_labels":["AIRPORT AREAS"],
                                "edge_categories_linewidth":None,
                                "edge_categories_zorder":None,
                                "node_categories_markersize":[15.0],
                                "node_categories_marker":["o"],
                                "node_categories_zorder":[11],
                            },
                            {
                                "sector":"transport",
                                "sector_gpkg":"", # We will call this separately
                                "sector_label":"Ports and Airports",
                                "edge_layer":None,
                                "node_layer":None,
                                "area_layer":None,
                                "edge_id_column":None,
                                "node_id_column":"node_id",
                                "area_id_column":None,
                                "edge_damage_filter_column":None,
                                "node_damage_filter_column":None,
                                "area_damage_filter_column":None,
                                "edge_damage_categories":None,
                                "node_damage_categories":None,
                                "area_damage_categories":None,
                                "edge_classify_column":None,
                                "node_classify_column":"asset_type",
                                "area_classify_column":None,
                                "edge_categories":None,
                                "node_categories":["port","airport"],
                                "area_categories":None,
                                "edge_categories_colors":None,
                                "node_categories_colors":["#a50f15","#8c510a"],
                                "area_categories_colors":None,
                                "edge_categories_labels":None,
                                "node_categories_labels":["PORT AREAS","AIRPORT AREAS"],
                                "area_categories_labels":None,
                                "edge_categories_linewidth":None,
                                "edge_categories_zorder":None,
                                "node_categories_markersize":[15.0,15.0],
                                "node_categories_marker":["o","o"],
                                "node_categories_zorder":[11,11],
                            },
                            {
                                "sector":"transport",
                                "sector_gpkg":"roads.gpkg",
                                "sector_label":"Roads",
                                "edge_layer":"edges",
                                "node_layer":"nodes",
                                "area_layer":None,
                                "edge_id_column":"edge_id",
                                "node_id_column":"node_id",
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
                                "edge_categories_colors":["#000000","#6a51a3","#ce1256","#f16913","#fdae6b","#fdae6b"],
                                "node_categories_colors":["#800026"],
                                "area_categories_colors":None,
                                "edge_categories_labels":["CLASS A","CLASS B","CLASS C","METRO","TRACK/OTHER","TRACK/OTHER"],
                                "node_categories_labels":["BRIDGES"],
                                "area_categories_labels":None,
                                "edge_categories_linewidth":[1.0,1.0,1.0,0.5,0.5,0.5],
                                "edge_categories_zorder":[10,9,8,7,6,6],
                                "node_categories_markersize":[10.0],
                                "node_categories_marker":".",
                                "node_categories_zorder":[11],
                                "node_loss_result": "road_bridges",
                                "edge_loss_result":"road_rail_edges"
                            },
                            # {
                            #     "sector":"energy",
                            #     "sector_gpkg":"electricity_network_v1.0.gpkg",
                            #     "sector_label":"Energy",
                            #     "edge_layer":"edges",
                            #     "node_layer":"nodes",
                            #     "area_layer":None,
                            #     "edge_id_column":"id",
                            #     "node_id_column":"id",
                            #     "area_id_column":None,
                            #     "edge_classify_column":"voltage",
                            #     "node_classify_column":"subtype",
                            #     "area_classify_column":None,
                            #     "edge_damage_filter_column":None,
                            #     "node_damage_filter_column":None,
                            #     "area_damage_filter_column":None,
                            #     "edge_damage_categories":None,
                            #     "node_damage_categories":None,
                            #     "area_damage_categories":None,
                            #     "edge_categories":["138kV","69 kV","24 kV","12 kV"],
                            #     "node_categories":["diesel","gas","hydro","solar","wind",
                            #                     "substation","pole","demand"],
                            #     "area_categories":None,
                            #     "edge_categories_colors":["#000000","#1f78b4","#33a02c","#f16913"],
                            #     "node_categories_colors":["#800026","#800026","#800026","#800026","#800026",
                            #                                 "#6a51a3","#ce1256","#f16913"],
                            #     "area_categories_colors":None,
                            #     "edge_categories_labels":["138kV","69kV","24kV","12kV"],
                            #     "node_categories_labels":["Power plant","Power plant",
                            #                             "Power plant","Power plant",
                            #                             "Power plant",
                            #                             "Substation",
                            #                             "Pole","Demand point"],
                            #     "area_categories_labels":None,
                            #     "edge_categories_linewidth":[1.0,1.0,1.0,1.0],
                            #     "edge_categories_zorder":[10,9,8,7],
                            #     "node_categories_markersize":[15.0,15.0,15.0,15.0,15.0,10.0,4.0,4.0],
                            #     "node_categories_marker":["s","s","s","s","s","o",".","."],
                            #     "node_categories_zorder":[20,20,20,20,20,19,18,17],
                            # },
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
                                # "node_categories":["diesel","gas","hydro","solar","wind",
                                #                 "substation","pole"],
                                "node_categories":["diesel","gas","hydro","solar","wind",
                                                "substation"],
                                "area_categories":None,
                                "edge_categories_colors":["#000000","#006d2c","#41ab5d","#66c2a4"],
                                # "node_categories_colors":["#800026","#800026","#800026","#800026","#800026",
                                #                             "#6a51a3","#ce1256"],
                                "node_categories_colors":["#800026","#800026","#800026","#800026","#800026",
                                                            "#6a51a3"],
                                "area_categories_colors":None,
                                "edge_categories_labels":["138kV lines","69kV lines","24kV lines","12kV lines"],
                                # "node_categories_labels":["Power plants","Power plants",
                                #                         "Power plants","Power plants",
                                #                         "Power plants",
                                #                         "Substations", 
                                #                         "Poles"
                                #                         ],
                                "node_categories_labels":["Power plants","Power plants",
                                                        "Power plants","Power plants",
                                                        "Power plants",
                                                        "Substations"
                                                        ],
                                "area_categories_labels":None,
                                "edge_categories_linewidth":[1.5,1.0,0.8,0.8],
                                "edge_categories_zorder":[10,9,8,7],
                                # "node_categories_markersize":[15.0,15.0,15.0,15.0,15.0,15.0,15.0],
                                # "node_categories_marker":["s","s","s","s","s","o","."],
                                # "node_categories_zorder":[15,15,15,15,15,14,13],
                                "node_categories_markersize":[15.0,15.0,15.0,15.0,15.0,15.0],
                                "node_categories_marker":["s","s","s","s","s","o"],
                                "node_categories_zorder":[15,15,15,15,15,14],
                                "node_loss_result": "electricity_nodes",
                                "edge_loss_result":"electricity_edges"
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
                                "edge_categories_colors":["#ec7014"],
                                "node_categories_colors":["#807dba",
                                                            "#6a51a3",
                                                            "#54278f",
                                                            "#3f007d",
                                                            "#4292c6",
                                                            "#2171b5",
                                                            "#08519c",
                                                            "#08306b",
                                                            "#41ab5d",
                                                            "#238b45",
                                                            "#006d2c",
                                                            "#00441b",
                                                            "#d94801",
                                                            "#a63603",
                                                            "#7f2704"
                                                            ],
                                "area_categories_colors":None,
                                "edge_categories_labels":["Pipelines"],
                                "node_categories_labels":["Booster Station","Catchment","Entombment",
                                                            "Filter Plant","Intake",
                                                            "Production Well","Pump Station",
                                                            "Relift Station","Reservoir",
                                                            "River Source","Spring","Storage Tank",
                                                            "Sump","Treatment Plant"],
                                "area_categories_labels":None,
                                "edge_categories_linewidth":[0.5],
                                "edge_categories_zorder":[10],
                                "node_categories_markersize":[6.0,6.0,6.0,
                                                                6.0,6.0,6.0,
                                                                6.0,6.0,6.0,
                                                                6.0,6.0,6.0,
                                                                6.0,6.0,6.0],
                                "node_categories_marker":["o","s","8",
                                                        "p","h","H",
                                                        "D","d","s",
                                                        ".","d","H",
                                                        "s","h","o"],
                                "node_categories_zorder":[15,15,15,
                                                            15,15,15,
                                                            15,15,15,
                                                            15,15,15,
                                                            15,15,15],
                                "node_loss_result": "potable_facilities",
                                "edge_loss_result":"potable_pipelines"                           
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
                                "edge_categories_colors":["#016c59","#3690c0"],
                                "node_categories_colors":["#014636"],
                                "area_categories_colors":None,
                                "edge_categories_labels":["Canals","Pipelines"],
                                "node_categories_labels":["Wells"
                                                        ],
                                "area_categories_labels":None,
                                "edge_categories_linewidth":[0.8,0.8],
                                "edge_categories_zorder":[10,10],
                                "node_categories_markersize":[10.0],
                                "node_categories_marker":["o"],
                                "node_categories_zorder":[15],
                                "node_loss_result": "irrigation_nodes",
                                "edge_loss_result":"irrigation_edges"
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
                                "edge_categories_colors":None,
                                "node_categories_colors":["#bfd3e6",
                                                    "#4d004b",
                                                    "#88419d",
                                                    "#ef3b2c"],
                                "area_categories_colors":None,
                                "edge_categories_labels":None,
                                "node_categories_labels":["Sump","Pumping Station","Relift Station","Treatment Plant"],
                                "area_categories_labels":None,
                                "edge_categories_linewidth":None,
                                "edge_categories_zorder":None,
                                "node_categories_markersize":[12.0,12.0,12.0,12.0],
                                "node_categories_marker":[".","^","s","o"],
                                "node_categories_zorder":[15,15,15,15],
                            },

                        ]

    return sector_attributes

def jamaica_port_and_airport_nodes():
    config = load_config()
    processed_data_path = config['paths']['data']


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

def jamaica_port_nodes():
    config = load_config()
    processed_data_path = config['paths']['data']

    transport_data_path = os.path.join(processed_data_path,
                        "networks",
                        "transport")
    ports = gpd.read_file(os.path.join(transport_data_path,"port_polygon.gpkg"),layer="areas").to_crs(epsg=JAMAICA_GRID_EPSG)
    ports["asset_type"] = "port"
    ports["geometry"] = ports.apply(lambda x: x.geometry.centroid,axis=1)
    ports = gpd.GeoDataFrame(ports,geometry="geometry",crs=f"EPSG:{JAMAICA_GRID_EPSG}")

    return ports
