"""Take the buildings footprints from OSM and add attributes to them
    Write a final buildings footprints into a geopackage
"""
import sys
import os

import geopandas as gpd
import pandas as pd
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
    vals = list(set([(x[cl[0]],x[cl[1]]) for cl in column_pairs]))
    if len(vals) > 1:
        vals = [v for v in vals if v[0] != 'X']
    
    return vals[0][0],vals[0][1]

def match_points_to_nearest_buildings(buildings_input,points_file,sector_codes):
    """Match points to their nearest buildings corresponding to a set of sector codes
    """
    buildings_select = buildings_input[buildings_input['sector_code'].isin(sector_codes)].reset_index()
    sindex_buildings = buildings_select.sindex
    points_file['bid'] = points_file.geometry.progress_apply(
        lambda x: get_nearest_node(x, sindex_buildings,buildings_select, 'bid'))
    points_file['geom'] = points_file.geometry.progress_apply(
        lambda x: get_nearest_node(x, sindex_buildings,buildings_select, 'geometry'))
    points_file['distance'] = points_file.progress_apply(lambda x: x.geometry.distance(x.geom),axis=1)
    points_file.drop(['geom'],axis=1,inplace=True)
    
    return points_file

def main(config):
    # Set the file paths on your machine for the data
    # All input data is stored in a path ../incoming_data/..
    # All output data will be written in a path ../processed_data/.. 
    
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    building_data_path = os.path.join(incoming_data_path,'buildings')
    epsg_jamaica = 3448  # Jamaica projectin system
    
    """Step 1: Extract the differet types of entries in useful columns 
        To infer what kind of building attributes are listed in the OSM building layer
        Using these attributes we want to map a macroeconomic sector to a building based on this attribute 
        The columns in the layer which have useful information are - 'building','amenity','office','shop','other_tags'
    """
    columns_types = ['building','amenity','office','shop','other_tags']
    buildings_input = gpd.read_file(os.path.join(building_data_path,'jamaica-buildings.gpkg'))
    buildings_input = buildings_input.to_crs(epsg=epsg_jamaica)
    buildings_input['bid'] = buildings_input.index.values.tolist()
    # print (buildings_input)
    """Extract categories within columns
    """
    output_excel = os.path.join(incoming_data_path,'buildings','unique_column_mapping.xlsx')
    output_wrtr = pd.ExcelWriter(output_excel) 
    for column in columns_types:
        df = buildings_input[column].value_counts().reset_index()
        df.columns = [column,'count']
        df.to_excel(output_wrtr,sheet_name=column, index=False)
        output_wrtr.save()
   

    """Step 2: Use a preprocessed file that maps buildings to economic sectors and subsectors
       The preprocessed file builds on the file output of Step 1 where now we have added two new columns
       sector_code & subsector_code - which are the codes for the macroeconomic sectors and subsectors respectively   
       The output of this step is to map the codes onto the spatial building layer to create two new columns in the layer 

       The default values for sector_code and subsector_code is 'X' - which signifies we cannot assign an economic attributes to the building 
    """
    merge_columns = []
    new_columns = []
    for column in columns_types:
        sector_map = pd.read_excel(os.path.join(
                                    building_data_path,
                                    'osm_building_column_mapping_economic_sectors.xlsx'),
                                    sheet_name=column)
        sector_map = sector_map[[column,'sector_code','subsector_code']]
        sector_map.rename(columns = {'sector_code':f'{column}_sector_code',
                                    'subsector_code':f'{column}_subsector_code'},
                        inplace=True)
        merge_columns += [(f'{column}_sector_code',f'{column}_subsector_code')]
        new_columns += [f'{column}_sector_code',f'{column}_subsector_code']
        buildings_input = pd.merge(buildings_input,sector_map,how='left',on=[column])
    
    buildings_input[new_columns] = buildings_input[new_columns].fillna(value='X')

    buildings_input['sector_subsector'] = buildings_input.progress_apply(lambda x: get_sector_subsector(x,merge_columns),axis=1)
    buildings_input[['sector_code','subsector_code']] = buildings_input['sector_subsector'].apply(pd.Series)
    buildings_input.drop(new_columns + ['sector_subsector'],axis=1,inplace=True)

    write_output = True
    if write_output is True: # We can write the intermediary output for sense checking if we want to
        buildings_input.to_file(os.path.join(
                                    building_data_path,
                                    'buildings_assigned_economic_sectors_intermediate.gpkg'), driver='GPKG')


    # print (buildings_input)
    # print (buildings_input.crs)
    
    """Step 3: Match buildings to known point locations provided from sources within Jamaica
        The known point locations are also mapped to macroeconomic sectors

        We find the nearest building to each point location 
        If that building is with a given distance buffer then we assign that point to the building

        We only match a point to a building which 
            either already has the same sector code assigned to it or is unassigned with code 'X' 
    """
    buildings_input['assigned_point'] = 'Unknown'
    distance_buffer = 1000 # buffer distance 
    buildings_input['centroid'] = buildings_input.progress_apply(lambda x:x.geometry.centroid,axis=1)
    # buildings_input.drop(['geometry'],axis=1,inplace=True)
    buildings_input.rename(columns={'geometry':'polygon_geometry'},inplace=True)
    buildings_input.rename(columns={'centroid':'geometry'},inplace=True)
    buildings_input = gpd.GeoDataFrame(buildings_input,geometry='geometry',crs={'init': f'epsg:{epsg_jamaica}'})

    points_layers_data = os.path.join(processed_data_path,'locations','points_of_interest.gpkg')

    nsdmb_layers = pd.read_excel(os.path.join(
                                    building_data_path,
                                    'nsdmb_layers_economic_sectors.xlsx'),
                                    sheet_name='Sheet1')
    for layer in nsdmb_layers.itertuples():
        layer_name = layer.layer
        sector_code = layer.sector_code
        subsector_code = layer.subsector_code
        points_file = gpd.read_file(points_layers_data,layer=layer_name)
        # print (points_file)
        # print (points_file.crs)
        points_file = points_file.to_crs(epsg=epsg_jamaica)

        points_file = match_points_to_nearest_buildings(buildings_input,points_file,['X',sector_code])
        points_file.to_csv(os.path.join(building_data_path,'test_values.csv'),index=False)
        building_ids = list(set(points_file[points_file['distance'] < distance_buffer]['bid'].values.tolist()))
        
        # print (building_ids)
        buildings_input.loc[buildings_input['bid'].isin(building_ids), 'sector_code'] = sector_code
        buildings_input.loc[buildings_input['bid'].isin(building_ids), 'subsector_code'] = subsector_code
        buildings_input.loc[buildings_input['bid'].isin(building_ids), 'assigned_point'] = layer_name

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
        points_file = gpd.read_file(points_layers_data,layer=layer['layer_name'])
        points_file = points_file.to_crs(epsg=epsg_jamaica)

        # This is to just fix a data problem in the indsutry layer
        if layer['layer_name'] == 'industry':
            points_file['SECTOR'].fillna('Unknown',inplace=True)
            for sh in range(len(layer['excel_sheets'])):
                nsdmb_layers = pd.read_excel(os.path.join(
                                    building_data_path,
                                    layer['excel_file']),
                                    sheet_name=layer['excel_sheets'][sh])
                points_file = pd.merge(points_file,nsdmb_layers,how='left',on=layer['layer_columns'][sh])
                points_file.rename(columns={'sector_code':f"{layer['excel_sheets'][sh]}_sector_code",
                                    'subsector_code':f"{layer['excel_sheets'][sh]}_subsector_code"},
                            inplace=True)

                points_file[f"{layer['excel_sheets'][sh]}_sector_code"].fillna('',inplace=True)
                points_file[f"{layer['excel_sheets'][sh]}_subsector_code"].fillna('',inplace=True)

            points_file['sector_code'] = points_file['manufacturing_sector_code'].astype(str) \
                                        + points_file['mining_sector_code'].astype(str)
            points_file['subsector_code'] = points_file['manufacturing_subsector_code'].astype(str) \
                                        + points_file['mining_subsector_code'].astype(str)
            points_file.drop(['manufacturing_sector_code','mining_sector_code',
                                'manufacturing_subsector_code','mining_subsector_code'],axis=1,inplace=True)
        else:
            nsdmb_layers = pd.read_excel(os.path.join(
                                    building_data_path,
                                    layer['excel_file'])
                                    )
            points_file = pd.merge(points_file,nsdmb_layers,how='left',on=layer['layer_columns'])

        points_sector_codes = ['X'] + list(set(points_file['sector_code'].values.tolist()))
        points_file = match_points_to_nearest_buildings(buildings_input,points_file,points_sector_codes)

        points_file = points_file[points_file['distance'] < distance_buffer]
        for point in points_file.itertuples():
            buildings_input.loc[buildings_input['bid'] == point.bid, 'sector_code'] = point.sector_code
            buildings_input.loc[buildings_input['bid'] == point.bid, 'subsector_code'] = point.subsector_code
            buildings_input.loc[buildings_input['bid'] == point.bid, 'assigned_point'] = layer['layer_name']

        print (f"* Done with points layer - {layer['layer_name']}")


    # Convert buuldings layer back to polygons and write to a geopackage
    buildings_input.drop(['geometry'],axis=1,inplace=True)
    buildings_input.rename(columns={'polygon_geometry':'geometry'},inplace=True)
    buildings_input = gpd.GeoDataFrame(buildings_input,geometry='geometry',crs={'init': f'epsg:{epsg_jamaica}'})

    write_output = True
    if write_output is True: # We can write the intermediary output for sense checking if we want to
        buildings_input.to_file(os.path.join(
                                    building_data_path,
                                    'buildings_assigned_economic_sectors.gpkg'), driver='GPKG')

if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)