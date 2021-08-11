"""Take different datasets from the NSDMB database and create geopackages
    Organise these geopackages into different types 
"""
import sys
import os

import geopandas as gpd
from preprocess_utils import *


def main(config):
    database_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]
    database_name = "GWP_Jamaica_NSP_Master_Geodatabase_v01.gdb"

    layer_details = [
        {
            "layer_type": "admin_boundaries",
            "layer_folder": "boundaries",
            "nsdmb_names": [
                "Administrative_Parishes",
                "Administrative_Communities_STATIN",
                "SociEconomic_Census_2011",
            ],
            "gpkg_names": ["admin1", "admin2", "admin3"],
        },
        {
            "layer_type": "land_type_and_use",
            "layer_folder": "land_type_and_use",
            "nsdmb_names": ["LandUse_LandUse", "soils", "Geology_Erosion_WMU"],
            "gpkg_names": ["land_use", "soils", "geology_erosion"],
        },
        {
            "layer_type": "port",
            "layer_folder": "networks",
            "nsdmb_names": ["shipping_ports"],
            "gpkg_names": ["ports"],
        },
        {
            "layer_type": "airport",
            "layer_folder": "networks",
            "nsdmb_names": ["airports_aerodromes_airstrips"],
            "gpkg_names": ["airports"],
        },
        {
            "layer_type": "points_of_interest",
            "layer_folder": "locations",
            "nsdmb_names": [
                "Community_Names",
                "community_centres",
                "health_centres",
                "industry",
                "markets",
                "Trade_Plant_Locations",
                "exporters",
                "bauxite_plants",
                "quarries",
                "postal_offices_agencies",
                "restaurants",
                "retail_stores",
                "cinemas",
                "court_houses",
                "CriticalFacilities_Cemetaries",
                "CriticalFacilities_Daycare",
                "fire_stations",
                "CriticalFacilities_PoliceStations",
                "CriticalFacilities_Schools",
                "libraries",
                "beaches_hotel_merge",
                "CriticalFacilities_Hotels",
                "hotels_guesthouses_villas",
                "cambios",
                "automotive_services",
                "car_dealers",
                "car_rental_agencies",
                "craft_markets",
                "CriticalFacilities_TouristAttractions",
                "food_packaging_houses",
            ],
            "gpkg_names": [
                "towns_and_villages",
                "community_centres",
                "health_centres",
                "industry",
                "markets",
                "trade_plants",
                "exporters",
                "bauxite_plants",
                "quarries",
                "postal_offices",
                "restaurants",
                "retail_stores",
                "cinemas",
                "court_houses",
                "cemetaries",
                "daycare",
                "firestations",
                "policestations",
                "schools",
                "libraries",
                "beach_hotels",
                "hotels",
                "hotels_guesthouses_villas",
                "money_exchange",
                "automotive_services",
                "car_dealers",
                "car_rental_agencies",
                "craft_markets",
                "touristattractions",
                "food_packaging_houses",
            ],
        },
        {
            "layer_type": "building_and_economic_activity_footprints",
            "layer_folder": "locations",
            "nsdmb_names": [
                "Hotosm_Jam_buildings_polygon",
                "aqua_farms",
                "fish_sanctuaries",
            ],
            "gpkg_names": ["buildings", "aqua_farms", "fish_sanctuaries"],
        },
    ]

    for layers in layer_details:
        output_path = os.path.join(processed_data_path, layers["layer_folder"])
        if os.path.exists(output_path) == False:
            os.mkdir(output_path)
        out_fname = os.path.join(output_path, f"{layers['layer_type']}.gpkg")
        for d in range(len(layers["nsdmb_names"])):
            nsdmb_layer = gpd.read_file(
                os.path.join(database_path, "nsdmb", database_name),
                layer=layers["nsdmb_names"][d],
            )
            nsdmb_layer.to_file(out_fname, layer=layers["gpkg_names"][d], driver="GPKG")
            print(f"* Done with layer {layers['nsdmb_names'][d]}")


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
