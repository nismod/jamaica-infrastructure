"""Create a alnd use planning database for Jamaica by extracting data from the NSDMB database 
    Map land use sector and subsector codes for the different land use types
"""
import os

import pandas as pd
import geopandas as gpd
import fiona
from collections import OrderedDict
from shapely.geometry import shape, mapping 
from preprocess_utils import *
from tqdm import tqdm
tqdm.pandas()

def get_sector_subsector_infra(x,column,class_dataframe):
    """This function filters out the most likely sector and subsector codes assigned to a building 
        In an OSM dataset, based on a mapping of building attributes to macroeconomic sectors

        Based on our mapping a building might be either assigned a sector based on some attributes 
            Or not assigned any sector, and instead assigned a code 'X' 
        If a building is assigned a code 'X' in addition to a known macroeconomic sector code
            Then the code 'X' is removed and the sector code is retained

        The output of the function results in 1 code assigned to a building   
    """
    landuse_split = pd.DataFrame(str(x[column]).split("/"),columns=["land_classes"])
    landuse_split = pd.merge(landuse_split,class_dataframe,how="left",on=["land_classes"])

    vals = list(OrderedDict.fromkeys(zip(landuse_split["sector_code"],
                                        landuse_split["subsector_code"])))
    if len(vals) > 1:
        sector_val = ",".join([str(v[0]) for v in vals if str(v[0]) != 'X'])
        subsector_val = ",".join([str(v[1]) for v in vals if str(v[1]) != 'X'])
    else:
        sector_val = vals[0][0]
        subsector_val = vals[0][1]
    
    vals = list(set(landuse_split["infra_sector"].values.tolist()))
    if len(vals) > 1:
        infra_val = ",".join([str(v) for v in vals if str(v) != 'X'])
    else:
        infra_val = vals[0]

    return sector_val, subsector_val, infra_val

def match_parishes_to_landplanning(jamaica_parishes,gdf,gdf_list):
    parish_match = gpd.sjoin(jamaica_parishes[["CODE","PARISH","geometry"]],
                                gdf[['layer_id','geometry']], 
                                how="inner", op='intersects').reset_index()
    # print (parish_match)
    parish_match.rename(columns={"geometry":"parish_geometry"},inplace=True)
    parish_match = pd.merge(parish_match, gdf[['layer_id','geometry']],how="left",on=["layer_id"])
    parish_match["area_match"] = parish_match.progress_apply(lambda x:x["parish_geometry"].intersection(x["geometry"].buffer(0)).area,axis=1)
    parish_match = parish_match.sort_values(by=["area_match"],ascending=False)
    parish_match = parish_match.drop_duplicates(subset=["layer_id"], keep="first")

    gdf = pd.merge(gdf,parish_match[["layer_id","PARISH"]],how="left",on=["layer_id"])
    parish = gdf["PARISH"].mode().values[0]
    gdf["PARISH"] = gdf["PARISH"].astype(str)
    gdf["PARISH"] = gdf["PARISH"].replace([r"^\s+$",'nan','None','none'], parish, regex=True)
    gdf.drop("layer_id",axis=1,inplace=True)
    gdf_list.append(gdf)

    return gdf_list

def main(config):
    incoming_data_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]
    epsg_jamaica = 3448
    
    """Step 1: Read the layers from the NSDMB database and write them into a Geopackage file
        This is a pre-preprocess step, so is commented out once it is done
    """
    # planning_layers = pd.read_excel(
    #                         os.path.join(
    #                             incoming_data_path,
    #                             "buildings",
    #                             "nsdmb_planning_layers.xlsx"),
    #                         sheet_name="Sheet1")
    # planning_layers["Type of Dataset"] = planning_layers.progress_apply(lambda x:x["Type of Dataset"].strip(),axis=1)
    # planning_layers = planning_layers[planning_layers["Type of Dataset"] == "Polygon"]

    # planning_layers.to_excel(os.path.join(
    #                             incoming_data_path,
    #                             "buildings",
    #                             "nsdmb_planning_layers_polygons_all.xlsx"),
    #                         sheet_name="Sheet1")
    
    # for p in planning_layers["Name of feature class "].values.tolist():
    #     layer_dict = []
    #     layer = fiona.open(os.path.join(
    #                             incoming_data_path,
    #                             "nsdmb",
    #                             "GWP_Jamaica_NSP_Master_Geodatabase_v01.gdb"), 
    #                         driver='OpenFileGDB', layer=p)
        
    #     for feature in layer:
    #         if any(feature["geometry"]["coordinates"]):
    #             d1 = {"geometry":shape(feature["geometry"])}
    #             d1.update(feature["properties"])
    #             layer_dict.append(d1)

    #     # Build the GeoDataFrame from Fiona Collection
    #     gdf = gpd.GeoDataFrame(pd.DataFrame(layer_dict),geometry="geometry", crs=layer.crs)
    #     gdf["layer_id"] = gdf.index.values.tolist()
    #     gdf.to_file(os.path.join(
    #                         incoming_data_path,
    #                         "buildings",
    #                         "planning_layers.gpkg"),layer=p, driver='GPKG')
    #     print ("* Done with layer",p)

    """Step 2: Take the excel file with layer names and query them from the database
    """
    planning_layers = pd.read_excel(os.path.join(
                                incoming_data_path,
                                "buildings",
                                "nsdmb_planning_layers_polygons.xlsx"),
                            sheet_name="Sheet1")
    columns = [c.strip() for c in planning_layers.columns.values.tolist()]
    planning_layers.columns = columns
    
    selected_layers = []
    # Clarendon land use layers
    gdf_merge = []
    for i, pl in planning_layers.iterrows():
        if "clarendon" in pl["Description"].lower() and "macarry" not in pl["Name of feature class"].lower():
            gdf = gpd.read_file(os.path.join(
                            incoming_data_path,
                            "buildings",
                            "planning_layers.gpkg"),layer=pl["Name of feature class"]).to_crs(epsg=epsg_jamaica)
            selected_layers.append(pl["Name of feature class"])
            if "LU_Zone" in gdf.columns.values.tolist():
                """find the most frequent land use and code and assign it to the blanks
                """
                for lu in ["LU_Zone","Code_Use"]:
                    lu_val = gdf[lu].mode().values[0]
                    gdf[lu] = gdf[lu].replace([r"^\s+$",'','nan','None','none'], lu_val, regex=True)

                gdf_merge.append(gdf[["LU_Zone","Code_Use","geometry"]])

        print ("* Done with layer",pl["Name of feature class"])

    gdf_merge = pd.concat(gdf_merge,axis=0,ignore_index=True)
    gdf_merge["land_use_id"] = gdf_merge.index.values.tolist()
    gdf = gpd.GeoDataFrame(gdf_merge,geometry="geometry", crs=f"EPSG:{epsg_jamaica}")
    gdf["area_sqm"] = gdf.progress_apply(lambda x:x.geometry.area,axis=1)
    gdf.to_file(os.path.join(
                        incoming_data_path,
                        "buildings",
                        "landuse_planning_layers.gpkg"),layer="clarendon_landuse", driver='GPKG')
    
    # gdf.groupby(["LU_Zone","Code_Use"])["area_sqm"].sum().reset_index().to_csv(os.path.join(
    #                             incoming_data_path,
    #                             "buildings",
    #                             "clarendon_landuse_planning_classes.csv"),index=False)

    # Manchester land use layers
    gdf_merge = []
    manchester_layers = pd.read_csv(os.path.join(
                            incoming_data_path,
                            "buildings","land_use_layers_uses_codes.csv"))
    for i, pl in planning_layers.iterrows():
        if "manchester" in pl["Description"].lower(): 
            ignore_layers = ["Kingston_Excavation_and_Reclaimation_Works",
                            "Portland_Blight_2007_Major_Land_Use",
                            "Goat_Island_Excavation_and_Reclaimation_Works",
                            "Entire_Portlandbight_PC_Division_1",
                            "Bowden_Excavation_and_Reclaimation_Works",
                            "land_use_Portland_proposed",
                            "Macarry_Bay_Excavation_and_Reclaimation_Works",
                            "Macarry_Bay_Option_2A_Excavation_and_Reclaimation_Works",
                            "NF_Macarry_Bay_Logistic_Hub_option_2B_detailed",
                            "NF_Macarry_Bay_option_2A_logistic_Hub_Concept_detailed",
                            "Water_Bodies_1"]
            if pl["Name of feature class"] not in ignore_layers:
                gdf = gpd.read_file(os.path.join(
                            incoming_data_path,
                            "buildings",
                            "planning_layers.gpkg"),layer=pl["Name of feature class"]).to_crs(epsg=epsg_jamaica)
                if pl["Name of feature class"] != "commercial_site_Manchester_proposed":
                    lu_zone = manchester_layers.loc[manchester_layers["layer_name"]== pl["Name of feature class"],"LU_Zone"].values[0]
                    lu_code = manchester_layers.loc[manchester_layers["layer_name"]== pl["Name of feature class"],"Code_Use"].values[0]
                    gdf["LU_Zone"] = lu_zone
                    gdf["Code_Use"] = lu_code
                else:
                    gdf["LU_Zone"] = "Commercial"
                
                gdf_merge.append(gdf[["LU_Zone","Code_Use","geometry"]])
                selected_layers.append(pl["Name of feature class"])

        print ("* Done with layer",pl["Name of feature class"])

    gdf_merge = pd.concat(gdf_merge,axis=0,ignore_index=True)
    gdf_merge["land_use_id"] = gdf_merge.index.values.tolist()
    gdf = gpd.GeoDataFrame(gdf_merge,geometry="geometry", crs=f"EPSG:{epsg_jamaica}")
    gdf["area_sqm"] = gdf.progress_apply(lambda x:x.geometry.area,axis=1)
    gdf.to_file(os.path.join(
                        incoming_data_path,
                        "buildings",
                        "landuse_planning_layers.gpkg"),
                        layer="manchester_landuse", 
                        driver='GPKG')
    # gdf.groupby(["LU_Zone","Code_Use"])["area_sqm"].sum().reset_index().to_csv(os.path.join(
    #                             incoming_data_path,
    #                             "buildings",
    #                             "manchester_landuse_planning_classes.csv"),index=False)
     
    # Create the landuse layers from the existing and proposed developments in Jamaica
    jamaica_parishes = gpd.read_file(os.path.join(
                                    processed_data_path,
                                    "boundaries",
                                    "admin_boundaries.gpkg"),
                            layer="admin1").to_crs(epsg=epsg_jamaica)
    existing_landuse_layers = ["ClarendonExistingLanduse",
                                "HanoverExistingLanduse",
                                "Negril_Existing_LandUse",
                                "PortmoreExistingLandUse",
                                "St_CatherineExistingLanduse",
                                "St_ElizabethExistingLanduse",
                                "St_MaryExistingLanduse",
                                "St_ThomasExistingLanduse",
                                "WestmorelandExistingLandUse",
                                "agriculture_SantaCruz",
                                "All_Jamaica_Feb_2018"]
    gdf_existing = []
    for layer in existing_landuse_layers:
        gdf = gpd.read_file(os.path.join(
                            incoming_data_path,
                            "buildings",
                            "planning_layers.gpkg"),layer=layer).to_crs(epsg=epsg_jamaica)
        gdf.rename(columns={"Existing_L":"LU_Zone",
                            "EX_LU":"LU_Zone",
                            "EX_LAND_US":"LU_Zone",
                            "LANDUSE":"LU_Zone",
                            "PROP_USE":"LU_Zone",},
                    inplace=True)
        gdf["LU_Zone"] = gdf["LU_Zone"].apply(str)
        gdf["LU_Zone"] = gdf.apply(lambda x: str(x["LU_Zone"]).strip(),axis=1)
        gdf["LU_Zone"] = gdf["LU_Zone"].replace([r"^\s+$",'','nan','None','none'], "Unknown", regex=True)
        if layer == "All_Jamaica_Feb_2018":
            gdf = gdf[gdf["LU_Zone"] != "Unknown"]
        gdf_existing = match_parishes_to_landplanning(jamaica_parishes,gdf[["layer_id","LU_Zone","geometry"]],gdf_existing)

        print ("* Done with layer",layer)
    
    gdf_existing = pd.concat(gdf_existing,axis=0,ignore_index=True)
    gdf_existing["land_use_id"] = gdf_existing.index.values.tolist()
    gdf = gpd.GeoDataFrame(gdf_existing,geometry="geometry", crs=f"EPSG:{epsg_jamaica}")
    gdf["area_sqm"] = gdf.progress_apply(lambda x:x.geometry.area,axis=1)
    gdf.to_file(os.path.join(
                        incoming_data_path,
                        "buildings",
                        "landuse_planning_layers.gpkg"),layer="existing_landuse", driver='GPKG')
    existing_classes = gdf.groupby(["LU_Zone"])["area_sqm"].sum().reset_index()

    proposed_landuse_layers = ["ClarendonProposals",
                                "HanoverProposals",
                                "land_use_Portland_proposed",
                                "PortmoreProposals",
                                "St_AnnExistingLandUse",
                                "St_CatherineProposals",
                                "St_ElizabethProposals",
                                "St_JamesProposals",
                                "St_MaryProposals",
                                "St_ThomasProposals",
                                "WestmorelandProposals",
                                "KSA_Proposals_2017_04_19"]
    gdf_existing = []
    for layer in proposed_landuse_layers:
        gdf = gpd.read_file(os.path.join(
                            incoming_data_path,
                            "buildings",
                            "planning_layers.gpkg"),layer=layer).to_crs(epsg=epsg_jamaica)
        # print (layer,gdf.columns.values.tolist())
        gdf.rename(columns={"Existing_L":"LU_Zone",
                            "EX_LU":"LU_Zone",
                            "EX_LAND_US":"LU_Zone",
                            "LANDUSE":"LU_Zone",
                            "PROP_USE":"LU_Zone",
                            "Pro_Use":"LU_Zone",
                            "pro_use":"LU_Zone"},
                    inplace=True)
        if "LU_Zone" not in gdf.columns.values.tolist():
            gdf["LU_Zone"] = 'Unknown'
        else:
            gdf["LU_Zone"] = gdf["LU_Zone"].apply(str)
            gdf["LU_Zone"] = gdf.apply(lambda x: str(x["LU_Zone"]).strip(),axis=1)
            gdf["LU_Zone"] = gdf["LU_Zone"].replace([r"^\s+$",'','','nan','None','none'], "Unknown", regex=True)
        
        if layer in ["land_use_Portland_proposed","St_AnnExistingLandUse"]:
            gdf.rename(columns={"Proposed_L":"Proposals"
                            },
                    inplace=True)

        if "Proposals" not in gdf.columns.values.tolist():
            gdf["Proposals"] = 'Unknown'
        else:
            gdf["Proposals"] = gdf["Proposals"].apply(str)
            gdf["Proposals"] = gdf.apply(lambda x: str(x["Proposals"]).strip(),axis=1)
            gdf["Proposals"] = gdf["Proposals"].replace([r"^\s+$",'','nan','None','none'], "Unknown", regex=True)
        
        gdf_existing = match_parishes_to_landplanning(jamaica_parishes,
                                    gdf[["layer_id","LU_Zone","Proposals","geometry"]],
                                    gdf_existing)
        
        print ("* Done with layer",layer)

    gdf_existing = pd.concat(gdf_existing,axis=0,ignore_index=True)
    gdf_existing["land_use_id"] = gdf_existing.index.values.tolist()
    gdf = gpd.GeoDataFrame(gdf_existing,geometry="geometry", crs=f"EPSG:{epsg_jamaica}")
    gdf["area_sqm"] = gdf.progress_apply(lambda x:x.geometry.area,axis=1)
    gdf.to_file(os.path.join(
                        incoming_data_path,
                        "buildings",
                        "landuse_planning_layers.gpkg"),layer="existing_proposed_landuse", driver='GPKG')

    # existing_classes = pd.concat([existing_classes, 
    #                             gdf.groupby(["LU_Zone"])["area_sqm"].sum().reset_index()],
    #                             axis=0,
    #                             ignore_index=True) 
    # existing_classes.groupby(["LU_Zone"])["area_sqm"].sum().reset_index().to_csv(os.path.join(
    #                             incoming_data_path,
    #                             "buildings",
    #                             "landuse_planning_existing_classes.csv"),index=False)

    # gdf.groupby(["Proposals"])["area_sqm"].sum().reset_index().to_csv(os.path.join(
    #                             incoming_data_path,
    #                             "buildings",
    #                             "landuse_planning_proposed_classes.csv"),index=False)

    # lu_classes = []
    # class_groups = pd.read_csv(os.path.join(
    #                             incoming_data_path,
    #                             "buildings",
    #                             "landuse_planning_existing_classes.csv"))
    # for i,cl in class_groups.iterrows():
    #     lu_classes += str(cl["LU_Zone"]).split("/")

    # class_groups = pd.read_csv(os.path.join(
    #                             incoming_data_path,
    #                             "buildings",
    #                             "landuse_planning_proposed_classes.csv"))
    # for i,cl in class_groups.iterrows():
    #     lu_classes += str(cl["Proposals"]).split("/")

    # lu_classes = list(set(lu_classes))
    # pd.DataFrame(lu_classes,columns=["land_classes"]).sort_values(by="land_classes").to_csv(os.path.join(
    #                             incoming_data_path,
    #                             "buildings",
    #                             "landuse_classes.csv"),index=False)

    """Step 3 map the classes to sector codes, subsector codes and infrastructure sectors
    """
    landuse_layers = ["clarendon_landuse",
                    "manchester_landuse",
                    "existing_landuse",
                    "existing_proposed_landuse"]
    landuse_classes = ["clarendon_landuse_planning_classes_with_sector_codes.csv",
                       "manchester_landuse_planning_classes_with_sector_codes.csv",
                       "landuse_classes_with_sector_codes.csv",
                       "landuse_classes_with_sector_codes.csv"]

    for i,(layer,sector_classes) in enumerate(list(zip(landuse_layers,landuse_classes))):
        gdf = gpd.read_file(os.path.join(
                                incoming_data_path,
                                "buildings",
                                "landuse_planning_layers.gpkg"),
                            layer=layer)
        csv = pd.read_csv(os.path.join(
                            incoming_data_path,
                            "buildings",
                            sector_classes))

        if layer in ["clarendon_landuse","manchester_landuse"]:
            csv.drop("area_sqm",axis=1,inplace=True)
            gdf = pd.merge(gdf,csv,how="left",on=["LU_Zone","Code_Use"])
        else:
            for col in ["LU_Zone","Proposals"]:
                if col in gdf.columns.values.tolist():
                    gdf["sector_subsector_infra"] =  gdf.progress_apply(
                                                        lambda x: get_sector_subsector_infra(x,
                                                                                    col,
                                                                                    csv),
                                                        axis=1)
                    if col == "LU_Zone":
                        gdf[["sector_code_existing",
                            "subsector_code_existing",
                            "infra_sector_existing"]] = gdf["sector_subsector_infra"].apply(pd.Series)
                    else:
                        gdf[["sector_code_proposals",
                            "subsector_code_proposals",
                            "infra_sector_proposals"]] = gdf["sector_subsector_infra"].apply(pd.Series)

                    gdf.drop("sector_subsector_infra",axis=1,inplace=True)

        gdf = gpd.GeoDataFrame(gdf,geometry="geometry", crs=f"EPSG:{epsg_jamaica}")
        gdf.to_file(os.path.join(
                    incoming_data_path,
                    "buildings",
                    "landuse_planning_layers_with_sectors.gpkg"),layer=layer, driver='GPKG')
        print ("* Done with layer",layer)

    
    # Remove the file landuse_planning_layers.gpkg because it is superseded by landuse_planning_layers_with_sectors.gpkg
    try:
        os.remove(os.path.join(
                    incoming_data_path,
                    "buildings",
                    "landuse_planning_layers.gpkg"))
    except OSError:
        pass 

if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
