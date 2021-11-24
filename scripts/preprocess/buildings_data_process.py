"""Take the buildings footprints from OSM and add attributes to them
    Write a final buildings footprints into a geopackage
"""
import sys
import os

import geopandas as gpd
import pandas as pd
import numpy as np
from collections import OrderedDict
from preprocess_utils import *
from tqdm import tqdm
tqdm.pandas()

def get_sector_subsector(x,column_pairs):
    """This function filters out the most likely sector and subsector codes assigned to a building 
        In an OSM dataset, based on a mapping of building attributes to macroeconomic sectors

        Based on our mapping a building might be either assigned a sector based on some attributes 
            Or not assigned any sector, and instead assigned a code 'X' 
        If a building is assigned a code 'X' in addition to a known macroeconomic sector code
            Then the code 'X' is removed and the sector code is retained

        The output of the function results in 1 code assigned to a building   
    """
    p0 = []
    p1 = []
    for cl in column_pairs:
        p0 += x[cl[0]].split(',')
        p1 += x[cl[1]].split(',')

    vals = list(OrderedDict.fromkeys(zip(p0,p1)))
    if len(vals) > 1:
        sector_val = ",".join([v[0] for v in vals if v[0] != 'X'])
        subsector_val = ",".join([v[1] for v in vals if v[1] != 'X'])
    else:
        sector_val = vals[0][0]
        subsector_val = vals[0][1]
    
    return sector_val, subsector_val

def join_values(x,column_types):
    """This function filters out the most likely sector and subsector codes assigned to a building 
        In an OSM dataset, based on a mapping of building attributes to macroeconomic sectors

        Based on our mapping a building might be either assigned a sector based on some attributes 
            Or not assigned any sector, and instead assigned a code 'X' 
        If a building is assigned a code 'X' in addition to a known macroeconomic sector code
            Then the code 'X' is removed and the sector code is retained

        The output of the function results in 1 code assigned to a building   
    """
    p0 = []
    for cl in column_types:
        p0 += str(x[cl]).split(',')

    vals = list(OrderedDict.fromkeys(p0))
    
    return ",".join([v for v in vals if v not in ['X','nan','None','none']])

def join_values_in_column(x,column_merge):
    """This function filters out the most likely sector and subsector codes assigned to a building 
        In an OSM dataset, based on a mapping of building attributes to macroeconomic sectors

        Based on our mapping a building might be either assigned a sector based on some attributes 
            Or not assigned any sector, and instead assigned a code 'X' 
        If a building is assigned a code 'X' in addition to a known macroeconomic sector code
            Then the code 'X' is removed and the sector code is retained

        The output of the function results in 1 code assigned to a building   
    """
    p0 = str(x[column_merge]).split(',')
    vals = list(OrderedDict.fromkeys(p0))
    
    return ",".join([v for v in vals if v not in ['X','nan','None','none']])

def create_sector_subsector_attribute_columns(building_dataframe,building_id_column,merge_columns,new_columns,columns_types):
    # building_dataframe[new_columns] = building_dataframe[new_columns].fillna('X',inplace=True)
    building_dataframe[new_columns] = building_dataframe[new_columns].astype('str')
    building_dataframe[new_columns] = building_dataframe[new_columns].replace([r"^\s+$",'nan','None','none'], 'X', regex=True)

    building_dataframe['sector_subsector'] = building_dataframe.progress_apply(lambda x: get_sector_subsector(x,merge_columns),axis=1)
    building_dataframe[['sector_code','subsector_code']] = building_dataframe['sector_subsector'].apply(pd.Series)
    columns_types = [c for c in columns_types if c != 'other_tags']
    building_dataframe['assigned_attribute'] = building_dataframe.progress_apply(
                                            lambda x:join_values(x,columns_types),
                                            axis=1)
    building_dataframe.drop(new_columns + columns_types + ['sector_subsector'],axis=1,inplace=True)
    building_dataframe = building_dataframe.drop_duplicates(subset=[building_id_column],keep='first')

    return building_dataframe

def sector_subsector_asset_mapping(building_dataframe,building_id_column,building_sector_mapping_file,columns_types):
    merge_columns = []
    new_columns = []
    for column in columns_types:
        sector_map = pd.read_excel(
                                    building_sector_mapping_file,
                                    sheet_name=column)
        sector_map = sector_map.drop_duplicates([column,'sector_code','subsector_code'],keep='first')
        sector_map = sector_map[[column,'sector_code','subsector_code']] 
        sector_map.rename(columns = {'sector_code':f'{column}_sector_code',
                                    'subsector_code':f'{column}_subsector_code'},
                        inplace=True)
        merge_columns += [(f'{column}_sector_code',f'{column}_subsector_code')]
        new_columns += [f'{column}_sector_code',f'{column}_subsector_code']
        building_dataframe = pd.merge(building_dataframe,sector_map,how='left',on=[column])

    return create_sector_subsector_attribute_columns(building_dataframe,building_id_column,merge_columns,new_columns,columns_types)

def match_points_to_nearest_buildings(points_of_interest,points_file,sector_codes):
    """Match points to their nearest buildings corresponding to a set of sector codes
    """
    buildings_select = points_of_interest[points_of_interest['sector_code'].isin(sector_codes)].reset_index()
    sindex_buildings = buildings_select.sindex
    points_file['osm_id'] = points_file.geometry.progress_apply(
        lambda x: get_nearest_node(x, sindex_buildings,buildings_select, 'osm_id'))
    points_file['geom'] = points_file.geometry.progress_apply(
        lambda x: get_nearest_node(x, sindex_buildings,buildings_select, 'geometry'))
    points_file['distance'] = points_file.progress_apply(lambda x: x.geometry.distance(x.geom),axis=1)
    points_file.drop(['geom'],axis=1,inplace=True)
    
    return points_file

def identify_bauxite_areas(x):
    if 'bauxite' in str(x["Classify"]).lower() or 'bauxite' in str(x["NAME"]).lower() or 'bauxite' in str(x["global_LU_type"]).lower():
        return 1
    else:
        return 0

def merge_columns(value_dataframe,building_dataframe,building_id_column,dataframe_merge_columns):
    for column in ['sector_code','subsector_code']:
        value_dataframe.rename(columns={column:f"{column}_match"},inplace=True)
        building_dataframe.rename(columns={column:f"{column}_building"},inplace=True)

    merge_columns = []
    new_columns = []
    # columns_types = []
    for column in ['match','building']:
        merge_columns += [(f'sector_code_{column}',f'subsector_code_{column}')]
        new_columns += [f'sector_code_{column}',f'subsector_code_{column}']
        # columns_types.append(f'assigned_attribute_{column}')

    building_dataframe = pd.merge(building_dataframe,value_dataframe,how='left',on=dataframe_merge_columns)
    building_dataframe[new_columns] = building_dataframe[new_columns].astype('str')
    building_dataframe[new_columns] = building_dataframe[new_columns].replace([r"^\s+$",'nan','None','none'], 'X', regex=True)
    building_dataframe['sector_subsector'] = building_dataframe.progress_apply(lambda x: get_sector_subsector(x,merge_columns),axis=1)
    building_dataframe[['sector_code','subsector_code']] = building_dataframe['sector_subsector'].apply(pd.Series)
    building_dataframe.drop(new_columns + ['sector_subsector'],axis=1,inplace=True)
    building_dataframe = building_dataframe.drop_duplicates(subset=[building_id_column],keep='first')

    return building_dataframe

def main(config):
    # Set the file paths on your machine for the data
    # All input data is stored in a path ../incoming_data/..
    # All output data will be written in a path ../processed_data/.. 
    
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    building_data_path = os.path.join(incoming_data_path,'buildings')
    points_of_interest_data_path = os.path.join(incoming_data_path,
                                        'hotosm_data',
                                        'hotosm_jam_points_of_interest_gpkg')
    epsg_jamaica = 3448  # Jamaica projection system
    
    # """Step 1: Extract the differet types of entries in useful columns 
    #     To infer what kind of building attributes are listed in the OSM building layer
    #     Using these attributes we want to map a macroeconomic sector to a building based on this attribute 
    #     The columns in the layer which have useful information are - 'building','amenity','office','shop','other_tags'
    # """
    # buildings_input = gpd.read_file(os.path.join(building_data_path,'jamaica-buildings.gpkg'))
    # buildings_input = buildings_input.to_crs(epsg=epsg_jamaica)
    # # buildings_input['osm_id'] = buildings_input.progress_apply(lambda x: modify_osm_id(x),axis=1)
    # buildings_input[['osm_id','osm_way_id']] = buildings_input[['osm_id','osm_way_id']].astype('str')
    # buildings_input[['osm_id','osm_way_id']] = buildings_input[['osm_id','osm_way_id']].replace([r"^\s+$",'nan','None','none'], '', regex=True)
    # buildings_input['osm_id'] = buildings_input['osm_id'].astype('str') + buildings_input['osm_way_id'].astype('str')

    # write_output = False
    # if write_output is True: # We can write the intermediary output for sense checking if we want to
    #     buildings_input.to_file(os.path.join(
    #                                 building_data_path,
    #                                 'buildings_assigned_economic_sectors_intermediate.gpkg'), driver='GPKG')


    # buildings_input_mapping = os.path.join(building_data_path,
    #                                         'osm_building_column_mapping_economic_sectors.xlsx')
    
    # points_of_interest = gpd.read_file(os.path.join(points_of_interest_data_path,
    #                                     'hotosm_jam_points_of_interest.gpkg'))
    # points_of_interest = points_of_interest.to_crs(epsg=epsg_jamaica)
    # points_of_interest['shop'] = points_of_interest['shop'].replace('yes','wholesale',regex=True)
    # points_of_interest_mapping = os.path.join(points_of_interest_data_path,
    #                                     'hotosm_mapping_economic_sectors.xlsx')

    # """
    # Excel outputs generated for post-processing and 
    # Creating a mapping of known building tags in the datasets with economic sectors
    # The rest of the code uses the outputs files from these tags


    # # Extract categories within columns
    
    # output_excel = os.path.join(incoming_data_path,'buildings','unique_column_mapping.xlsx')
    # output_wrtr = pd.ExcelWriter(output_excel) 
    # for column in columns_types:
    #     df = buildings_input[column].value_counts().reset_index()
    #     df.columns = [column,'count']
    #     df.to_excel(output_wrtr,sheet_name=column, index=False)

    # output_wrtr.save()

    # # Extract categories within columns
    
    # columns_types = ['building','amenity','office','shop']
    # sector_maps = []
    # for column in columns_types:
    #     sector_map = pd.read_excel(os.path.join(
    #                                 building_data_path,
    #                                 'osm_building_column_mapping_economic_sectors.xlsx'),
    #                                 sheet_name=column)
    #     sector_map = sector_map[[column,'sector_code','subsector_code']]
    #     sector_map.rename(columns = {column:'asset_type'},inplace=True)
    #     sector_maps.append(sector_map)

    # sector_maps = pd.concat(sector_maps,axis=0,ignore_index=True)


    # columns_types = ['amenity','man_made','shop','tourism']
    # output_excel = os.path.join(points_of_interest_data_path,
    #                                 'unique_column_mapping.xlsx')
    # output_wrtr = pd.ExcelWriter(output_excel) 
    # for column in columns_types:
    #     df = points_of_interest[column].value_counts().reset_index()
    #     df.columns = [column,'count']
    #     df = pd.merge(df,sector_maps,how='left',left_on=[column],right_on='asset_type')
    #     df.drop('asset_type',axis=1,inplace=True)
    #     df.to_excel(output_wrtr,sheet_name=column, index=False)
    
    # output_wrtr.save()
    # """

    # """Step 2: Use a preprocessed file that maps buildings to economic sectors and subsectors
    #    The preprocessed file builds on the file output of Step 1 where now we have added two new columns
    #    sector_code & subsector_code - which are the codes for the macroeconomic sectors and subsectors respectively   
    #    The output of this step is to map the codes onto the spatial building layer to create two new columns in the layer 

    #    The default values for sector_code and subsector_code is 'X' - which signifies we cannot assign an economic attributes to the building 
    # """

    # points_of_interest = sector_subsector_asset_mapping(points_of_interest,
    #                                             'osm_id',
    #                                             points_of_interest_mapping,
    #                                             ['amenity','man_made','shop','tourism'])
    # write_output = False
    # if write_output is True: # We can write the intermediary output for sense checking if we want to
    #     points_of_interest.to_file(os.path.join(
    #                                 points_of_interest_data_path,
    #                                 'points_of_interest_assigned_economic_sectors_intermediate.gpkg'), driver='GPKG')
    # buildings_input = sector_subsector_asset_mapping(buildings_input,
    #                                             'osm_id',
    #                                             buildings_input_mapping,
    #                                             ['building','amenity','office','shop','other_tags'])
    # write_output = False
    # if write_output is True: # We can write the intermediary output for sense checking if we want to
    #     buildings_input.to_file(os.path.join(
    #                                 building_data_path,
    #                                 'buildings_assigned_economic_sectors_intermediate.gpkg'), driver='GPKG')

    
    # points_of_interest = points_of_interest.to_crs(epsg=epsg_jamaica)
    # buildings_input = buildings_input.to_crs(epsg=epsg_jamaica)
    
    # if 'osm_id' in points_of_interest.columns.values.tolist():
    #     points_of_interest.drop("osm_id",axis=1,inplace=True)

    # points_intersections = gpd.sjoin(points_of_interest,
    #                                 buildings_input[['osm_id','geometry']], 
    #                                 how="inner", op='intersects').reset_index()
    # df = pd.DataFrame(list(set(points_intersections['osm_id'].values.tolist())),columns=['osm_id'])

    # for column in ['sector_code','subsector_code','assigned_attribute']:
    #     df = pd.merge(df,points_intersections[['osm_id',
    #                                 column,
    #                                 ]].groupby(['osm_id'])[column].apply(','.join).reset_index(),
    #                                 how='left',on=['osm_id'])
    #     df.rename(columns={column:f"{column}_match"},inplace=True)
    #     buildings_input.rename(columns={column:f"{column}_building"},inplace=True)

    # merge_columns = []
    # new_columns = []
    # columns_types = []
    # for column in ['match','building']:
    #     merge_columns += [(f'sector_code_{column}',f'subsector_code_{column}')]
    #     new_columns += [f'sector_code_{column}',f'subsector_code_{column}']
    #     columns_types.append(f'assigned_attribute_{column}')

    # buildings_input = pd.merge(buildings_input,df,how='left',on=['osm_id'])
    # buildings_input = create_sector_subsector_attribute_columns(buildings_input,'osm_id',
    #                                                     merge_columns,
    #                                                     new_columns,
    #                                                     columns_types)

    # write_output = False
    # if write_output is True: # We can write the intermediary output for sense checking if we want to
    #     buildings_input.to_file(os.path.join(
    #                                 building_data_path,
    #                                 'buildings_assigned_economic_sectors_intermediate.gpkg'), driver='GPKG')

    # """Identify the buildings within mining and quarrying areas
    # """
    # # buildings_input = gpd.read_file(os.path.join(
    # #                                 building_data_path,
    # #                                 'buildings_assigned_economic_sectors_intermediate.gpkg'))

    # land_use_types = gpd.read_file(os.path.join(
    #                             processed_data_path,
    #                             "land_type_and_use",
    #                             "jamaica_land_use_combined.gpkg"))
    # land_use_types = land_use_types.to_crs(epsg=epsg_jamaica)
    # buildings_input = buildings_input.to_crs(epsg=epsg_jamaica)

    # if "index_left" in land_use_types.columns.values.tolist():
    #     land_use_types.drop("index_left",axis=1,inplace=True)
    # elif "index_right" in land_use_types.columns.values.tolist():
    #     land_use_types.drop("index_right",axis=1,inplace=True)

    # land_use_types['bauxite_identify'] = land_use_types.progress_apply(lambda x: identify_bauxite_areas(x),axis=1)
    # mining_quarry_areas = land_use_types[land_use_types['bauxite_identify'] == 1]

    # mining_intersections = gpd.sjoin(mining_quarry_areas,
    #                                 buildings_input[['osm_id','geometry']], 
    #                                 how="inner", op='intersects').reset_index()
    # df = pd.DataFrame(list(set(mining_intersections['osm_id'].values.tolist())),columns=['osm_id'])
    # df["sector_code"] = "C"
    # df["subsector_code"] = 132
    # df["assigned_attribute"] = "Bauxite mining and alumina"

    # for column in ['sector_code','subsector_code','assigned_attribute']:
    #     df.rename(columns={column:f"{column}_match"},inplace=True)
    #     buildings_input.rename(columns={column:f"{column}_building"},inplace=True)

    # merge_columns = []
    # new_columns = []
    # columns_types = []
    # for column in ['match','building']:
    #     merge_columns += [(f'sector_code_{column}',f'subsector_code_{column}')]
    #     new_columns += [f'sector_code_{column}',f'subsector_code_{column}']
    #     columns_types.append(f'assigned_attribute_{column}')

    # buildings_input = pd.merge(buildings_input,df,how='left',on=['osm_id'])
    # buildings_input = create_sector_subsector_attribute_columns(buildings_input,'osm_id',
    #                                                     merge_columns,
    #                                                     new_columns,
    #                                                     columns_types)
    # write_output = True
    # if write_output is True: # We can write the intermediary output for sense checking if we want to
    #     buildings_input.to_file(os.path.join(
    #                                 building_data_path,
    #                                 'buildings_assigned_economic_sectors_intermediate.gpkg'), driver='GPKG')

    """Add the commercial buildings mapping
    """
    commercial_buildings = gpd.read_file(os.path.join(
    							incoming_data_path,
    							"JAM_classified_buildings_industries",
    							"JAM_classified_buildings.shp"))[["osm_id","layer",
                                                                "Name_of_Pl","SECTOR1",
                                                                "SUB_SECTOR","PRODUCTS",
                                                                "HS_CODE","COMPANY_NA","geometry"]]
    commercial_buildings.rename(columns={"SECTOR1":"SECTOR","Name_of_Pl":"Name_of_Plant"},inplace=True)
    commercial_buildings["layer"] = commercial_buildings["layer"].astype('str')
    commercial_buildings["layer"] = commercial_buildings["layer"].replace([r"^\s+$",'nan','None','none',"points_of_interest",' '],'', regex=True)
    # print (list(set(commercial_buildings["layer"].values.tolist())))
    commercial_buildings['sector_code'] = 'X'
    commercial_buildings['subsector_code'] = 'X'
    nsdmb_layers = pd.read_excel(os.path.join(
                                    building_data_path,
                                    'nsdmb_layers_economic_sectors.xlsx'),
                                    sheet_name='Sheet1')
    commercial_layers = []
    for layer in nsdmb_layers.itertuples():
        layer_name = layer.layer
        sector_code = layer.sector_code
        subsector_code = layer.subsector_code
        
        commercial_layer = commercial_buildings[commercial_buildings["layer"] == layer_name]
        commercial_layer["sector_code"] = sector_code
        commercial_layer["subsector_code"] = subsector_code
        commercial_layers.append(commercial_layer)
        # print (building_ids)
        # commercial_buildings.loc[commercial_buildings['layer'] == layer_name, 'sector_code'] = sector_code
        # commercial_buildings.loc[commercial_buildings['layer'] == layer_name, 'subsector_code'] = subsector_code
        # commercial_buildings.loc[buildings_input['layer'] == layer_name, 'assigned_point'] = layer_name

        print (f'* Done with points layer - {layer_name}')

    """Step 4: Do similar point layer mapping for three specific layers of 
        Trade plants, Exporters, and Industry layers
    """

    layer_dict = [
                    {
                        'layer_name':'trade_plants',
                        'excel_file':'nsdmb_tradeplants_layer_economic_sectors.xlsx',
                        'layer_columns':['Name_of_Plant'],
                    },
                    {
                        'layer_name':'exporters',
                        'excel_file':'nsdmb_exporters_layer_economic_sectors.xlsx',
                        'layer_columns':['SECTOR','SUB_SECTOR','PRODUCTS','HS_CODE'],

                    },
                    {
                        'layer_name':'industry',
                        'excel_file':'nsdmb_industry_layer_economic_sectors.xlsx',
                        'excel_sheets':['manufacturing','mining'],
                        'layer_columns':[['SECTOR'],['SECTOR','COMPANY_NA']]
                    },
                ]

    for layer in layer_dict:
        # points_file = gpd.read_file(points_layers_data,layer=layer['layer_name'])
        # points_file = points_file.to_crs(epsg=epsg_jamaica)

        # This is to just fix a data problem in the indsutry layer
        commercial_layer = commercial_buildings[commercial_buildings["layer"] == layer['layer_name']]
        if layer['layer_name'] == 'industry':
            # points_file['SECTOR'] = points_file['SECTOR'].replace(' ', 'Unknown')
            for sh in range(len(layer['excel_sheets'])):
                nsdmb_layers = pd.read_excel(os.path.join(
                                    building_data_path,
                                    layer['excel_file']),
                                    sheet_name=layer['excel_sheets'][sh])
                commercial_layer = merge_columns(nsdmb_layers,commercial_layer,"osm_id",layer['layer_columns'][sh])
        else:
            nsdmb_layers = pd.read_excel(os.path.join(
                                    building_data_path,
                                    layer['excel_file'])
                                    )
            commercial_layer = merge_columns(nsdmb_layers,commercial_layer,"osm_id",layer['layer_columns'])

        commercial_layers.append(commercial_layer)
        print (f"* Done with points layer - {layer['layer_name']}")

    commercial_layers = pd.concat(commercial_layers,axis=0,ignore_index=True)
    commercial_buildings = commercial_buildings[~commercial_buildings["osm_id"].isin(commercial_layers["osm_id"].values.tolist())]
    commercial_buildings = pd.concat([commercial_layers,commercial_buildings],axis=0,ignore_index=True)
    with_geom = commercial_buildings[["osm_id","geometry"]].drop_duplicates(subset=["osm_id"],keep="first")
    without_geom = commercial_buildings.drop("geometry",axis=1)
    print (without_geom)
    without_geom = without_geom.set_index('osm_id').astype(str).groupby('osm_id').agg(',' .join).reset_index()

    mcs = [c for c in without_geom.columns.values.tolist() if c != "osm_id"]
    for mc in mcs:
        without_geom[mc] = without_geom.progress_apply(lambda x:join_values_in_column(x,mc),axis=1)

    # without_geom = without_geom.groupby('osm_id').agg(lambda x: ','.join(set(x))).reset_index()
    commercial_buildings = gpd.GeoDataFrame(pd.merge(without_geom,with_geom,
                                    how="left",on=["osm_id"]),
                                    geometry="geometry",
                                    crs=f"EPSG:{epsg_jamaica}")
    write_output = True
    if write_output is True: # We can write the intermediary output for sense checking if we want to
        commercial_buildings.to_file(os.path.join(
                                    building_data_path,
                                    'commercial_buildings_assigned_economic_sectors_intermediate.gpkg'), driver='GPKG')

    # residential_buildiings = gpd.read_file(os.path.join(
    #                                 incoming_data_path,
    #                                 "JAM_residential",
    #                                 "jam_residential.shp"))


if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)