"""Create a road network for Jamaica
    Take road data from the NSDMB database and create topological networks
    Match the geometries and attributes between a NWA national roads network
    With more detailed road network
    Add the final network into a geopackage

"""
import sys
import os

import pandas as pd
import geopandas as gpd
import numpy as np
from shapely.geometry import Point,LineString
import assign_network_components
from collections import defaultdict
from preprocess_utils import *
from tqdm import tqdm
tqdm.pandas()

def modify_road_class(x):
    if x['first_clas'] in [None,np.nan]:
        x['first_clas'] = ''
    if x['road_class'] in [None,np.nan]:
        x['road_class'] = ''

    if x['first_clas']:
        return str(x['first_clas'])
    elif x['road_class']:
        return str(x['road_class'])
    else:
        return 'METRO'

def remodify_road_class(x):
    if x["road_class"] not in [None,np.nan,' ','']:
        return x["road_class"]
    elif x["tag_highway"] == "trunk":
        return "CLASS A"
    else:
    	return str(x["tag_highway"]).upper()

def modify_road_surface(x):
    if x['constructi'] in [None,np.nan,' ']:
        x['constructi'] = ''
    if x['construc_1'] in [None,np.nan,' ']:
        x['construc_1'] = ''

    x['constructi'] = str(x['constructi'])
    x['construc_1'] = str(x['construc_1'])

    if (x['constructi']) or (x['construc_1']):
        if x['constructi'] == x['construc_1']:
            return x['constructi']
        elif (x['constructi'] in ['SD & AC','AC & SD']) and (x['construc_1'] in ['SD & AC','AC & SD']):
            return 'SD & AC'
        elif x['construc_1']:
            return x['construc_1']
        else:
            return x['constructi']
    else:
        return 'Surface Dressed'

def remodify_road_surface(x):
    if x['road_surface'] in ['SD & AC','AC & SD']:
        return 'SD & AC'
    elif x['road_surface'] not in [None,np.nan,'',' ']:
    	return x["road_surface"]
    elif x["material"] in ['asphalt','as1']:
        return "Asphaltic Concrete" 
    else:
        return 'Surface Dressed'    

def modify_road_surface_again(x):
    if x['road_surface'] not in [None,np.nan,' ']:
        return str(x['road_surface'])
    else:
        return 'Surface Dressed' 

def modify_street_name(x):
    if x['street_name'] not in [None,np.nan,' ']:
        return str(x['street_name'])
    else:
        return x['str_name']

def modify_road_width(x,standard_width=6.1):
    if float(x['averagewid']) > 0:
        return float(x['averagewid'])
    elif float(x['average_wi']) > 0:
        return float(x['average_wi'])
    else:
        return standard_width

def remodify_road_width(x,standard_width=3.05):
    if float(x['road_width']) > 0:
        return float(x['road_width'])
    else:
        return x["lanes"]*standard_width

def modify_traffic_count(x):
    if float(x['aadt_1']) > 0:
        return float(x['aadt_1'])
    elif float(x['aadt']) > 0:
        return float(x['aadt'])
    else:
        return 0 

def modify_road_section(x):
    if x['section_1'] in [None,np.nan,' ']:
        x['section_1'] = ''
    if x['section'] in [None,np.nan,' ']:
        x['section'] = ''        
    if str(x['section_1']):
        return x['section_1']
    elif str(x['section']):
        return x['section']
    else:
        return 'unknown'

def remodify_road_section(x):
    if x['section_name'] in [None,np.nan,' ']:
    	return x['section_name']
    else:
        return 'unknown'

def modify_traffic_count_with_points(x):
    if x['Total_IN'] > 0:
        return float(x['Total_IN'])
    else:
        return x['traffic_count']

def modify_road_section_with_points(x):
    if x['NWA Section No.'] not in [None,np.nan,'']:
        return x['NWA Section No.']
    else:
        return x['section_name']

def modify_road_class_with_points(x):
    if x['NWA Road Class'] not in [None,np.nan,'']:
        return x['NWA Road Class']
    else:
        return x['road_class']

def modify_merged_sections(x):
    if x['NWA Section No.'] in [None,np.nan]:
        return ''
    else:
        return x['NWA Section No.']

def modify_merged_class(x):
    if x['NWA Road Class'] in [None,np.nan,'']:
        return ''
    else:
        return x['NWA Road Class']

def create_road_asset_type(x):
    if x["asset_type"] == "road_bridge":
    	return x["asset_type"]
    elif x["road_surface"] in ['Asphaltic Concrete', 'SD & AC', 'Surface Dressed']:
        return 'road_paved'
    else:
        return 'road_unpaved'

def match_roads(buffer_dataframe,edge_dataframe,
                geom_buffer=10,fraction_intersection=0.95,length_intersected=100,save_buffer_file=False):
    # print (nwa_edges)
    buffer_dataframe['geometry'] = buffer_dataframe.geometry.progress_apply(lambda x: x.buffer(geom_buffer))
    
    # Save the result to sense check by visual inspection on QGIS. Not a necessary step 
    if save_buffer_file is not False:
        buffer_dataframe.to_file(save_buffer_file,layer=f'nwa_buffer_{geom_buffer}',driver='GPKG')

    road_matches = gpd.sjoin(edge_dataframe,buffer_dataframe[['nwa_edge_id','geometry']], how="inner", op='intersects').reset_index()
    if len(road_matches.index) > 0:
        buffer_dataframe.rename(columns={'geometry':'nwa_geometry'},inplace=True)
        road_matches = pd.merge(road_matches,buffer_dataframe[['nwa_edge_id','nwa_geometry','nwa_length']],how='left',on=['nwa_edge_id'])
        
        # Find the length intersected and its percentage as the length of the road segment and the NWA road segment
        road_matches['length_intersected'] = road_matches.progress_apply(
                                lambda x: (x.geometry.intersection(x.nwa_geometry).length),
                                axis=1)
        road_matches['fraction_intersection'] = road_matches.progress_apply(
                                lambda x: (x.geometry.intersection(x.nwa_geometry).length)/x.geometry.length,
                                axis=1)
        road_matches['fraction_buffer'] = road_matches.progress_apply(
                                lambda x: (x.geometry.intersection(x.nwa_geometry).length)/x['nwa_length'],
                                axis=1)

        road_matches.drop(['nwa_geometry'],axis=1,inplace=True)
        
        # Filter out the roads whose 95%(0.95) or over 100-meters length intersects with the buffer  
        return road_matches[
                        (road_matches['fraction_intersection']>=fraction_intersection
                        ) | (road_matches['length_intersected']>=length_intersected)]
    else:
        return pd.DataFrame()

# Get the road network id's and the NWA road ID's that uniquely match
def get_unique_matches(dataframe):   
    df = dataframe['edge_id'].value_counts().reset_index()
    df.columns = ['edge_id','count']
    return df, dataframe[dataframe['edge_id'].isin(df[df['count'] == 1]['edge_id'].values.tolist())]

def road_length_matches_filtering(remaining_matches,
                        buffer_dataframe,edge_dataframe,
                        geom_buffer,geom_buffer_divisor,save_buffer_file=False):
    if len(remaining_matches.index) > 0:
        matches = []
        while len(remaining_matches.index) > 0:
            geom_buffer = 1.0*geom_buffer/geom_buffer_divisor
            print ('* Longest length matches remaining',len(remaining_matches.index))
            filter_remaining_matches = remaining_matches.groupby('edge_id')['fraction_intersection'].max().reset_index()
            filter_remaining_matches['longest_section_match'] = 1
            remaining_matches = pd.merge(remaining_matches,filter_remaining_matches,how='left',on=['edge_id','fraction_intersection']).fillna(0)
            remaining_matches = remaining_matches[remaining_matches['longest_section_match'] == 1]
            df1, unique_remaining_matches = get_unique_matches(remaining_matches)
            matches.append(unique_remaining_matches)

            print ('* Unique matches done',len(unique_remaining_matches))
            remaining_matches = remaining_matches[remaining_matches['edge_id'].isin(df1[df1['count'] > 1]['edge_id'].values.tolist())]
            print ('* Matches left',len(remaining_matches.index))
            if len(remaining_matches.index) > 0:
                nwa_ids = list(set(remaining_matches['nwa_edge_id'].values.tolist()))
                edge_ids = list(set(remaining_matches['edge_id'].values.tolist()))
                remaining_matches = match_roads(buffer_dataframe[buffer_dataframe['nwa_edge_id'].isin(nwa_ids)].copy(),
                                  edge_dataframe[edge_dataframe['edge_id'].isin(edge_ids)].copy(),
                                  geom_buffer=geom_buffer,
                                  fraction_intersection=0.5,
                                  save_buffer_file=save_buffer_file)
            # print ('* Longest length matches remaining',len(remaining_matches.index))
        remaining_matches = pd.concat(matches,axis=0,ignore_index=True)
    
    return remaining_matches

def assign_node_asset_type(x):
    if x['PARISH'] not in [None,np.nan,'']:
        return 'bridge'
    else:
        return 'junction'

def add_edge_speeds(x):
    if x.tag_maxspeed > 0:
    	return (float(x.tag_maxspeed))
    elif x.road_class in ["CLASS A","CLASS B","MOTORWAY","PRIMARY","SECONDARY"]:
        return 110
    elif x.road_class in ["CLASS C","TERTIARY"]:
        return 80
    else:
        return 50

def add_road_lanes(x):
    if x["road_class"] in ["CLASS A","CLASS B","CLASS C"]:
        lanes = round(x["road_width"]/3.65,0) 
    else:
        lanes = round(x["road_width"]/3.048,0)
    if lanes % 2 == 0:
        return lanes
    else:
        return lanes + 1

def remodify_road_lanes(x,standard_width=3.05):
    lanes = round(x["width_m"]/standard_width,0)
    if lanes % 2 == 0:
        return lanes
    else:
        return lanes + 1

def main(config):
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    epsg_jamaica = 3448

    road_network_path = os.path.join(incoming_data_path,'roads')
    save_intermediary_results = True # If you want to store the intermediary results for sense checking
    store_intersections = os.path.join(road_network_path,'road_intersections.gpkg')

    """
    Step 0: Get a processed and clean OSM dervied network layer
    """
    implement_step = False
    if implement_step is True:
	    edges = gpd.read_parquet(os.path.join(road_network_path,
	                                "different_roads",
	                                "edges.geoparquet"))    
	    nodes = gpd.read_parquet(os.path.join(road_network_path,
	                                "different_roads",
	                                "nodes.geoparquet"))
	    network = add_network_topology(nodes,edges,"road")
	    edges = network.edges
	    nodes = network.nodes
	    edges = edges.to_crs(epsg=epsg_jamaica)
	    nodes = nodes.to_crs(epsg=epsg_jamaica)

	    largest_component = edges['component_id'].value_counts().idxmax()
	    nodes = nodes[nodes.component_id == largest_component]
	    edges = edges[edges.component_id == largest_component]

	    if save_intermediary_results is True:
	        edges.to_file(os.path.join(road_network_path,"osm_roads.gpkg"),layer='edges',driver='GPKG')
	        nodes.to_file(os.path.join(road_network_path,"osm_roads.gpkg"),layer='nodes',driver='GPKG')
    
    """
    Step 1: Match two NWA layers and merge their properties to finalise an NWA layer
    """
    implement_step = True
    if implement_step is True:
	    nwa_edges = gpd.read_file(os.path.join(road_network_path,'different_roads','roads_main_NWA_NSDMD.gpkg'))
	    nwa_edges.columns = map(str.lower, nwa_edges.columns)
	    nwa_edges["averagewid"] = nwa_edges["averagewid"].fillna(0)
	    nwa_edges["average_wi"] = nwa_edges["average_wi"].replace(' ','0')
	    nwa_edges["road_surface"] = nwa_edges.progress_apply(lambda x:modify_road_surface(x),axis=1)
	    nwa_edges['section_name'] = nwa_edges.progress_apply(lambda x:modify_road_section(x),axis=1)
	    nwa_edges['traffic_count'] = nwa_edges.progress_apply(lambda x:modify_traffic_count(x),axis=1)
	    nwa_edges["average_wi"] = nwa_edges.progress_apply(lambda x:modify_road_width(x),axis=1)

	    nwa_edges = nwa_edges.drop_duplicates(subset=["section_name"],keep="first")
	    nwa_roads = gpd.read_file(os.path.join(road_network_path,'NWA','main_road_NWA.shp'))
	    nwa_roads = nwa_roads.to_crs(epsg=epsg_jamaica)

	    # Deal with empty edges (drop)
	    empty_idx = nwa_roads.geometry.apply(lambda e: e is None or e.is_empty)
	    if empty_idx.sum():
	        empty_edges = nwa_roads[empty_idx]
	        print(f"Found {len(empty_edges)} empty edges.")
	        print(empty_edges)
	        nwa_roads = nwa_roads[~empty_idx]
	        del empty_edges, empty_idx

	    nwa_roads.columns = map(str.lower, nwa_roads.columns)
	    nwa_roads.rename(columns={"section_":"section_name","street_nam":"street_name"},inplace=True)
	    nwa_roads = pd.merge(nwa_roads,
	                    nwa_edges[["section_name","str_name",
	                                "road_surface","average_wi",
	                                "traffic_count"]],how="left",on=["section_name"])
	    nwa_roads["street_name"] = nwa_roads.progress_apply(lambda x:modify_street_name(x),axis=1)
	    nwa_roads["road_width"] = nwa_roads.progress_apply(lambda x:modify_road_width(x),axis=1)
	    nwa_roads["road_surface"] = nwa_roads.progress_apply(lambda x:modify_road_surface_again(x),axis=1)
	    nwa_roads["traffic_count"] = nwa_roads["traffic_count"].fillna(0)

	    del nwa_edges
	    nwa_roads = nwa_roads[["objectid","section_name",
	                        "road_class","street_name",
	                        "road_surface","road_width","traffic_count","geometry"]]
	    # nwa_roads['nwa_length'] = nwa_roads.geometry.length
	    # nwa_roads.rename(columns={'objectid':'nwa_edge_id'},inplace=True)
	    # remove same roads, which might be creating issues
	    nwa_roads = nwa_roads[~nwa_roads["objectid"].isin([145,147,112,149])]
	    nwa_roads = nwa_roads.explode()
	    # nwa_roads["start"] = nwa_roads.progress_apply(lambda x: Point(x.geometry.coords[0]),axis=1)
	    # nwa_roads["end"] = nwa_roads.progress_apply(lambda x: Point(x.geometry.coords[-1]),axis=1)
	    nwa_roads = nwa_roads.reset_index()
	    nwa_roads["nwa_edge_id"] = nwa_roads.index.values.tolist()
	    nwa_roads["nwa_length"] = nwa_roads.geometry.length

    # Most of the geometries between the two networks are within a 5-meter buffer of each other
    # Match the two networks by creating a 20-meter buffer around the NWA roads and intersecting with the roads networks
    # We also select a road if it intersects more than 100-meters of the NWA buffer
    """Alternative solution
    """
    implement_step = False
    if implement_step is True:
	    edges = gpd.read_file(os.path.join(road_network_path,"osm_roads.gpkg"),layer='edges')
	    nodes = gpd.read_file(os.path.join(road_network_path,"osm_roads.gpkg"),layer='nodes')

	    nwa_edges = nwa_roads.copy()
	    road_select = match_roads(nwa_edges,edges,geom_buffer=5)
	    save_intermediary_results = False
	    if save_intermediary_results is True:
	        road_select.to_file(store_intersections,layer='selected_roads',driver='GPKG')

	    # Get the road network id's and the NWA road ID's that uniquely match    
	    df1, unique_matches = get_unique_matches(road_select)
	    if save_intermediary_results is True:
	        unique_matches.to_file(store_intersections,layer='unique_matches',driver='GPKG')

	    print ('* Total unique matches - first pass:',len(unique_matches.index))
	    print ('* Total unique ID matches - first pass:',len(list(set(unique_matches['edge_id'].values.tolist()))))
	    multiple_matches = road_select[road_select['edge_id'].isin(df1[df1['count'] > 1]['edge_id'].values.tolist())]
	    multiple_matches['length_difference'] = multiple_matches.progress_apply(lambda x:abs(x['fraction_intersection'] - x['fraction_buffer']),axis=1)
	    closest_matches = multiple_matches.groupby('edge_id')['length_difference'].min().reset_index()
	    closest_matches = closest_matches[closest_matches['length_difference'] <= 1e-2]
	    closest_matches['closeness'] = 1
	    multiple_matches = pd.merge(multiple_matches,closest_matches,how='left',on=['edge_id','length_difference']).fillna(0)
	    filter_multiple_matches = multiple_matches[multiple_matches['closeness'] == 1]
	    if save_intermediary_results is True:
	        filter_multiple_matches.to_file(store_intersections,layer='multiple_matches_filter',driver='GPKG')

	    print ('* Total equal length matches - first pass:',len(filter_multiple_matches.index))
	    print ('* Total equal length ID matches - first pass:',len(list(set(filter_multiple_matches['edge_id'].values.tolist()))))

	    finished_matches = list(set(unique_matches['edge_id'].values.tolist() + filter_multiple_matches['edge_id'].values.tolist()))
	    remaining_matches = road_select[~road_select['edge_id'].isin(finished_matches)]
	    remaining_matches = road_length_matches_filtering(remaining_matches,
	                        nwa_roads,edges,
	                        10,1.5,save_buffer_file=False)

	    if save_intermediary_results is True:
	        remaining_matches.to_file(store_intersections,
	                                                    layer='length_matches_filter',driver='GPKG')

	    print ('* Total longest length of intersection matches - first pass:',len(remaining_matches.index))
	    print ('* Total longest length of intersection matches ID matches - first pass:',
	                    len(list(set(remaining_matches['edge_id'].values.tolist()))))

	    all_matched_pairs = [unique_matches[['edge_id','nwa_edge_id']],
	                            filter_multiple_matches[['edge_id','nwa_edge_id']],
	                            remaining_matches[['edge_id','nwa_edge_id']]]
	    all_matches = list(set(unique_matches['edge_id'].values.tolist() \
	                    + filter_multiple_matches['edge_id'].values.tolist() \
	                    + remaining_matches['edge_id'].values.tolist()))
	    if len(road_select[~road_select['edge_id'].isin(all_matches)]) > 0:
	        print ('* Some roads still not matched')
	        print (road_select[~road_select['edge_id'].isin(all_matches)])
	    else:
	        print ('* First pass of matching done')

	    """Step 3: Check for the remaining unmatched NWA roads
	        Match them based on intersections with wider buffers 
	        And find the longest lengths of intersections
	    """
	    nwa_column_check = ['nwa_edge_id','fraction_buffer']
	    nwa_matches = pd.concat([unique_matches[nwa_column_check],
	                            filter_multiple_matches[nwa_column_check],
	                            remaining_matches[nwa_column_check]],
	                            axis=0,ignore_index=True)
	    nwa_matches = nwa_matches.groupby('nwa_edge_id')['fraction_buffer'].sum().reset_index()
	    nwa_edges = pd.merge(nwa_roads,nwa_matches,how='left',on=['nwa_edge_id']).fillna(0)
	    if save_intermediary_results is True:
	        nwa_edges[nwa_edges['fraction_buffer'] <= 0.5].to_file(store_intersections,layer='nwa_remaining_matches',driver='GPKG')

	    nwa_edges = nwa_edges[nwa_edges['fraction_buffer'] <= 0.5]
	    match_edges = edges[~edges['edge_id'].isin(all_matches)]
	    # match_edges = match_edges[~match_edges['_7'].isin(['OTHER','TRACK'])] # There are a lot of unwanted edges if we do not filter out this case
	    road_select = match_roads(nwa_edges,
	                        match_edges,
	                        geom_buffer=30,
	                        fraction_intersection=0.7,
	                        length_intersected=120,
	                        save_buffer_file=False)
	    if save_intermediary_results is True:
	        road_select.to_file(store_intersections,layer='wider_buffer_match',driver='GPKG')

	    remaining_matches = road_length_matches_filtering(road_select,
	                        nwa_roads,edges,
	                        30,1.01,save_buffer_file=False)
	    if save_intermediary_results is True:
	        remaining_matches.to_file(store_intersections,
	                                layer='length_matches_filter_wider_buffer',driver='GPKG')

	    print ('* Total longest length of intersection matches - second pass:',len(remaining_matches.index))
	    print ('* Total longest length of intersection matches ID matches - second pass:',
	                    len(list(set(remaining_matches['edge_id'].values.tolist()))))

	    all_matched_pairs += [remaining_matches[['edge_id','nwa_edge_id']]]

	    all_matched_pairs = pd.concat(all_matched_pairs,axis=0,ignore_index=True)

	    print ('* Total matches:',len(all_matched_pairs.index))
	    print ('* Total ID matches',len(list(set(all_matched_pairs['edge_id'].values.tolist()))))
	    print (all_matched_pairs)

	    save_intermediary_results = True
	    if save_intermediary_results is True:
	        df = gpd.GeoDataFrame(pd.merge(edges,all_matched_pairs,how='left',on=['edge_id']),
	                        geometry="geometry",crs=f"EPSG:{epsg_jamaica}")
	        df.to_file(store_intersections,layer='final_matches',driver='GPKG')
	        # nwa_roads[~nwa_roads['nwa_edge_id'].isin(
	        #             list(set(all_matched_pairs['nwa_edge_id'].values.tolist()))
	        #             )].to_file(store_intersections,layer='unmatched',driver='GPKG')
    
    """Step 4: Now match the attributes between the two networks
    """
    # Select the NWA attributes we want to retain
    nwa_columns = ["nwa_edge_id",
                    "section_name",
                    "road_class","street_name",
                    "road_surface","road_width",
                    "traffic_count" 
                    ]
    # all_matched_pairs = pd.merge(all_matched_pairs,nwa_roads[nwa_columns],how='left',on=['nwa_edge_id'])
    edges = gpd.read_file(store_intersections,layer='final_matches')
    edges = pd.merge(edges,nwa_roads[nwa_columns],how='left',on=['nwa_edge_id'])
    # Modify the edges columns first
    # edges = edges.rename(
    #     columns={
    #         "_7": "road_class",
    #         "street_nam": "street_name",
    #         "street_typ": "street_type",
    #         "section":'section_name'
    #     }
    # ).drop(
    #     columns=[
    #         "fnode_",
    #         "tnode_",
    #         "lpoly_",
    #         "rpoly_",
    #         "length",
    #         "shape_length",
    #         "road50west",
    #         "road50we_1",
    #     ]
    # )
    # edges = pd.merge(edges,all_matched_pairs,how='left',on=['edge_id'])
    # print (edges.columns)

    # string_columns = ['road_class','first_clas',
    #                 'street_name','str_name',
    #                 'section_name','section','section_1',
    #                 'constructi','construc_1',
    #                 'street_type','vertalignm']
    # edges[string_columns] = edges[string_columns].replace(r'^\s*$', '', regex=True)
    # # edges[string_columns] = edges[string_columns].fillna('',inplace=True)
    # numeric_columns = ['aadt','aadt_1','averagewid','average_wi']
    # for num in numeric_columns:
    #     edges[num] = edges[num].replace(' ', '0.0')
    #     edges[num] = edges[num].astype(float)
    #     # print ('* Converted column',num)

    # edges.to_file(os.path.join(road_network_path,'roads.gpkg'),layer='edges_modified_intermediate',driver="GPKG")

    # # edges = gpd.read_file(os.path.join(road_network_path,'roads.gpkg'),layer='edges_modified_intermediate')

    edges['road_class'] = edges.progress_apply(lambda x:remodify_road_class(x),axis=1)
    edges['road_surface'] = edges.progress_apply(lambda x:remodify_road_surface(x),axis=1)
    edges['section_name'] = edges.progress_apply(lambda x:remodify_road_section(x),axis=1)
    edges["lanes"] = np.where(edges["lanes"] == 1,2,edges["lanes"])
    # edges['street_name'] = edges.progress_apply(lambda x:modify_street_name(x),axis=1)
    edges['width_m'] = edges.progress_apply(lambda x:remodify_road_width(x),axis=1)
    edges["lanes"] = edges.progress_apply(lambda x:remodify_road_lanes(x),axis=1)
    # edges['traffic_count'] = edges.progress_apply(lambda x:modify_traffic_count(x),axis=1)

    # # print (edges.columns)
    # edges = edges.rename(
    #     columns={
    #         "from_id": "from_node",
    #         "to_id": "to_node",
    #         "id": "edge_id",
    #     }
    # ).drop(
    #     columns=[
    #         'nwa_edge_id',
    #         'first_clas',
    #         'str_name',
    #         'constructi',
    #         'averagewid', 
    #         'section', 
    #         'aadt',
    #         'section_1', 
    #         'average_wi', 
    #         'construc_1', 
    #         'aadt_1'
    #     ]
    # )
    edges = gpd.GeoDataFrame(edges,geometry="geometry",crs=f"EPSG:{epsg_jamaica}")
    edges.to_file(store_intersections,layer='edges_attributes',driver="GPKG")

    """Step 5: Modify the road traffic counts with the NWA points data
    """
    # edges = gpd.read_file(os.path.join(road_network_path,'roads.gpkg'),layer='edges_modified',driver="GPKG")
    edges.set_crs(epsg=epsg_jamaica)
    traffic_points = pd.read_csv(os.path.join(incoming_data_path,
                                    'roads','nwa_traffic_counts',
                                    'Updated NWA Traffic Count_July 27 2021.csv'))
    traffic_points = traffic_points[['NWA Section No.','NWA Road Class','Total_IN','X Coordinate', 'Y Coordinate']]
    traffic_points['point_id'] = traffic_points.index.values.tolist()
    geometry = [Point(xy) for xy in zip(traffic_points['X Coordinate'], traffic_points['Y Coordinate'])]
    traffic_points = gpd.GeoDataFrame(traffic_points, crs={'init': 'epsg:4326'}, geometry=geometry)
    traffic_points = traffic_points.to_crs(epsg=epsg_jamaica)

    # print (nwa_edges)
    geom_buffer = 100
    edges_buffered = edges.copy()
    edges_buffered['geometry'] = edges_buffered.geometry.progress_apply(lambda x: x.buffer(geom_buffer))
    point_matches = gpd.sjoin(traffic_points,edges_buffered[['edge_id','section_name','geometry']], how="inner", op='intersects').reset_index()

    point_matches.rename(columns = {'geometry':'point_geometry'},inplace=True)
    point_matches = pd.merge(point_matches,edges[['edge_id','geometry']],how='left',on=['edge_id'])
    point_matches['distance'] = point_matches.progress_apply(lambda x: x.geometry.distance(x.point_geometry),axis=1)
    # point_matches.to_csv(os.path.join(incoming_data_path,'roads','nwa_traffic_counts','matched_roads.csv'),index=False)

    all_matches = []
    same_matches = point_matches[point_matches['NWA Section No.'] == point_matches['section_name']]
    nearest_points = same_matches.groupby('point_id')['distance'].min().reset_index()
    nearest_points = pd.merge(nearest_points,
                same_matches[['point_id','edge_id','NWA Section No.','NWA Road Class','section_name','Total_IN','distance']],
                how='left',
                on=['point_id','distance'])
    all_matches.append(nearest_points)
    # point_ids = list(set(nearest_points['point_id'].values.tolist()))
    remaining_matches = point_matches[~point_matches['point_id'].isin(list(set(nearest_points['point_id'].values.tolist())))]
    nearest_points = remaining_matches.groupby('point_id')['distance'].min().reset_index()
    nearest_points = pd.merge(nearest_points,
                remaining_matches[['point_id','edge_id','NWA Section No.','NWA Road Class','section_name','Total_IN','distance']],
                how='left',
                on=['point_id','distance'])
    all_matches.append(nearest_points)
    all_matches = pd.concat(all_matches,axis=0,ignore_index=True)
    df = all_matches['edge_id'].value_counts().reset_index()
    df.columns = ['edge_id','count']
    all_matches = pd.merge(all_matches,df,how='left',on=['edge_id'])
    # all_matches.to_csv(os.path.join(incoming_data_path,'roads','nwa_traffic_counts','nearest_roads_final.csv'),index=False) 
    all_matches['NWA Section No.'] = all_matches.apply(lambda x:modify_merged_sections(x),axis=1)
    all_matches['NWA Road Class'] = all_matches.apply(lambda x:modify_merged_class(x),axis=1)

    all_matches_totals = all_matches.groupby('edge_id')['Total_IN'].mean().reset_index()
    all_matches_sections = all_matches.groupby('edge_id')['NWA Section No.'].apply(list).reset_index()
    all_matches_sections['NWA Section No.'] = all_matches_sections.apply(lambda x: ','.join(list(set(x['NWA Section No.']))),axis=1)
    all_matches_class = all_matches.groupby('edge_id')['NWA Road Class'].apply(list).reset_index()
    all_matches_class['NWA Road Class'] = all_matches_class.apply(lambda x: ''.join(list(set(x['NWA Road Class']))),axis=1)

    all_matches = pd.merge(all_matches_sections,all_matches_totals,how='left',on=['edge_id'])
    all_matches = pd.merge(all_matches,all_matches_class,how='left',on=['edge_id'])
    edges = pd.merge(edges,all_matches,how='left',on=['edge_id'])
    edges['section_name'] = edges.progress_apply(lambda x:modify_road_section_with_points(x),axis=1)
    edges['traffic_count'] = edges.progress_apply(lambda x:modify_traffic_count_with_points(x),axis=1)
    edges['traffic_count'] = edges['traffic_count'].fillna(0)
    edges['road_class'] = edges.progress_apply(lambda x:modify_road_class_with_points(x),axis=1)

    # print (edges)
    edges = edges.drop(columns=['NWA Section No.','NWA Road Class','Total_IN'])
    edges["tag_maxspeed"]  = edges["tag_maxspeed"].fillna(0.0)
    edges["speed_kph"] = edges.progress_apply(lambda x:add_edge_speeds(x),axis=1)
    edges["asset_type"] = edges.progress_apply(lambda x:create_road_asset_type(x),axis=1)
    # edges.to_file(os.path.join(processed_data_path,'networks','transport','roads.gpkg'),layer='edges',driver="GPKG")
    # edges.to_file(store_intersections,layer='edges_final_attributes',driver="GPKG")

    # nodes = gpd.read_file(os.path.join(road_network_path,'roads.gpkg'),layer='nodes')
    # nodes.set_crs(epsg=epsg_jamaica)
    # nodes['asset_type'] = nodes.progress_apply(lambda x:assign_node_asset_type(x),axis=1)
    # nodes.to_file(os.path.join(processed_data_path,'networks','transport','roads.gpkg'),layer='nodes',driver="GPKG")
    # From NWA - Cost of 2 lane road reconstruction is US$ 1.5 million/km
    # We convert it into US$ 0.75 million/km/lane and to J$ by assuming exchange rate is 1 J$ = 0.0068 US$
    road_costs = 0.5*1.5*1.0e6/1.0e3/0.0068
    bridge_costs = 1.5*1.0e6/0.0068
    # print (costs)

    edges = gpd.GeoDataFrame(edges,geometry="geometry",crs=f"EPSG:{epsg_jamaica}")
    edges["length_m"] = edges.geometry.length
    edges["cost_unit"] = "J$/m"
    edges["mean_damage_cost"] = np.where(edges["asset_type"] == "road_bridge",
    										bridge_costs/edges["length_m"],
    										road_costs*edges["lanes"])
    edges["min_damage_cost"] = 0.8*edges["mean_damage_cost"]
    edges["max_damage_cost"] = 1.2*edges["mean_damage_cost"]
    edges = gpd.GeoDataFrame(edges,geometry="geometry",crs=f"EPSG:{epsg_jamaica}")
    edges.to_file(store_intersections,layer='edges_final_attributes',driver="GPKG")
    # edges["reopen_cost"] = 152000.0 # Estimate from NWA
    # edges["reopen_cost_unit"] = "J$/day"
    # edges.drop(['min_reopen_cost',
    #         'mean_reopen_cost',
    #         'max_reopen_cost'],axis=1,inplace=True)
    edges.rename(columns={"road"})
    edges = edges[["edge_id","from_node","to_node",
    			"road_class","street_name",
    			"section_name","asset_type",
    			"road_surface","lanes",
    			"width_m","length_m","speed_kph",
    			"cost_unit","min_damage_cost",
    			"mean_damage_cost","max_damage_cost",
    			"component_id",
    			"geometry"]]
    edges.to_file(os.path.join(processed_data_path,
                                    'networks',
                                    'transport',
                                    'roads.gpkg'),
                            layer='edges')
    nodes = gpd.read_file(os.path.join(road_network_path,"osm_roads.gpkg"),layer='nodes')
    # nodes = gpd.read_file(os.path.join(processed_data_path,
    #                                 'networks',
    #                                 'transport',
    #                                 'roads.gpkg'),
    #                         layer='nodes')
    # nodes.rename(columns={"id":"node_id"},inplace=True)
    # nodes.drop(['min_reopen_cost',
    #         'mean_reopen_cost',
    #         'max_reopen_cost'],axis=1,inplace=True)
    # bridges = nodes[nodes["asset_type"] == "bridge"]
    # non_bridges = nodes[nodes["asset_type"] != "bridge"]
    # bridges["reopen_cost"] = 152000.0 # Estimate from NWA
    # bridges["reopen_cost_unit"] = "J$/day" 
    
    # costs = 1.5*1.0e6/0.0067
    # bridges["mean_damage_cost"] = costs
    # bridges["min_damage_cost"] = 0.8*costs
    # bridges["max_damage_cost"] = 1.2*costs

    # nodes = gpd.GeoDataFrame(pd.concat([bridges,non_bridges],axis=0,ignore_index=True),
    #             geometry="geometry",crs=f"EPSG:{epsg_jamaica}")
    # print (nodes)

    nodes.to_file(os.path.join(processed_data_path,
                                    'networks',
                                    'transport',
                                    'roads.gpkg'),
                            layer='nodes')
if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)