"""Create a road network for Jamaica
    Take road data from the NSDMB database and create topological networks
    Match the geometries and attributes between a NWA national rads network
    With more detailed road network
    Add the final network into a geopackage

"""
import sys
import os

import geopandas as gpd
from preprocess_utils import *
from tqdm import tqdm
tqdm.pandas()

def match_intersection_fractions(x):
    if abs(x['fraction_intersection'] - x['fraction_buffer']) <= 1.0e-2:
        return 1
    else:
        return 0

def match_roads(nwa_edges,edges,geom_buffer=10,fraction_intersection=0.95,length_intersected=100,save_buffer_file=False):
    nwa_edges['geometry'] = nwa_edges.geometry.progress_apply(lambda x: x.buffer(geom_buffer))
    
    # Save the result to sense check by visual inspection on QGIS. Not a necessary step 
    if save_buffer_file is not False:
        nwa_edges.to_file(save_buffer_file,layer=f'nwa_buffer_{geom_buffer}',driver='GPKG')

    road_matches = gpd.sjoin(edges,nwa_edges[['nwa_edge_id','geometry']], how="inner", op='intersects').reset_index()
    nwa_edges.rename(columns={'geometry':'nwa_geometry'},inplace=True)
    road_matches = pd.merge(road_matches,nwa_edges[['nwa_edge_id','nwa_geometry','nwa_length']],how='left',on=['nwa_edge_id'])
    
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
    road_select = road_matches[
                    (road_matches['fraction_intersection']>=fraction_intersection
                    ) | (road_matches['length_intersected']>=length_intersected)]
    return road_select

# Get the road network id's and the NWA road ID's that uniquely match
def get_unique_matches(road_select):   
    df = road_select['id'].value_counts().reset_index()
    df.columns = ['id','count']
    return df, road_select[road_select['id'].isin(df[df['count'] == 1]['id'].values.tolist())]

def main(config):
    database_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    epsg_jamaica = 3448

    """
    Step 1: Take the NWA layer and convert it into a tologogical network
    """
    # nwa_edges = gpd.read_file(os.path.join(database_path,'nsdmb','GWP_Jamaica_NSP_Master_Geodatabase_v01.gdb'),
    #                     layer='roads_main_NWA')
    # create_network_from_nodes_and_edges(None,nwa_edges,
    #                                 'nwa_road',
    #                                 os.path.join(processed_data_path,'networks','transport',"nwa_roads.gpkg"),by=None)

    """
    Step 2: Take another detailed network with preprocessed topology creation
    Match the road network to the NWA roads to get prperties of roads
    """
    edges = gpd.read_file(os.path.join(processed_data_path,'networks','transport','roads.gpkg'),layer='edges')
    edges.set_crs(epsg=epsg_jamaica)
    
    nwa_roads = gpd.read_file(os.path.join(processed_data_path,'networks','transport','nwa_roads.gpkg'),layer='edges')
    nwa_roads.set_crs(epsg=epsg_jamaica)
    nwa_roads['nwa_length'] = nwa_roads.progress_apply(lambda x: x['geometry'].length,axis=1)
    nwa_roads = nwa_roads[nwa_roads['nwa_length']>0]
    nwa_roads.rename(columns={'edge_id':'nwa_edge_id'},inplace=True)
    
    # Most of the geometries between the two networks are within a 10-meter buffer of each other
    # Match the two networks by creating a 10-meter buffer around the NWA roads and intersecting with the roads networks
    # We also select a road if it intersects mre than 100-meters of the NWA buffer
    
    store_intersections = os.path.join(processed_data_path,'networks','transport','road_intersections.gpkg')
    
    nwa_edges = nwa_roads.copy()
    road_select = match_roads(nwa_edges,edges)
    road_select.to_file(store_intersections,layer='selected_roads',driver='GPKG')

    # Get the road network id's and the NWA road ID's that uniquely match    
    df1, unique_matches = get_unique_matches(road_select)
    unique_matches.to_file(store_intersections,layer='unique_matches',driver='GPKG')

    multiple_matches = road_select[road_select['id'].isin(df1[df1['count'] > 1]['id'].values.tolist())]
    multiple_matches['length_difference'] = multiple_matches.progress_apply(lambda x:abs(x['fraction_intersection'] - x['fraction_buffer']),axis=1)
    closest_matches = multiple_matches.groupby('id')['length_difference'].min().reset_index()
    closest_matches = closest_matches[closest_matches['length_difference'] <= 1e-2]
    closest_matches['closeness'] = 1
    multiple_matches = pd.merge(multiple_matches,closest_matches,how='left',on=['id','length_difference']).fillna(0)
    filter_multiple_matches = multiple_matches[multiple_matches['closeness'] == 1]
    filter_multiple_matches.to_file(store_intersections,layer='multiple_matches_filter',driver='GPKG')

    finished_matches = unique_matches['id'].values.tolist() + filter_multiple_matches['id'].values.tolist()
    remaining_matches = road_select[~road_select['id'].isin(finished_matches)]
    if len(remaining_matches.index) > 0:
        filter_remaining_matches = remaining_matches.groupby('id')['fraction_intersection'].max().reset_index()
        filter_remaining_matches['longest_section_match'] = 1
        remaining_matches = pd.merge(remaining_matches,filter_remaining_matches,how='left',on=['id','fraction_intersection']).fillna(0)
        remaining_matches[remaining_matches['longest_section_match'] == 1].to_file(store_intersections,
                                                    layer='length_matches_filter',driver='GPKG')

    all_matches = list(set(unique_matches['id'].values.tolist() \
                    + filter_multiple_matches['id'].values.tolist() \
                    + remaining_matches['id'].values.tolist()))
    if len(road_select[~road_select['id'].isin(all_matches)]) > 0:
        print ('* Some roads still not matched')
        print (road_select[~road_select['id'].isin(all_matches)])
    else:
        print ('* First pass of matching done')

    # Check for thhe reamining unmatched NWA roads
    nwa_column_check = ['nwa_edge_id','fraction_buffer']
    nwa_matches = pd.concat([unique_matches[nwa_column_check],
                            filter_multiple_matches[nwa_column_check],
                            remaining_matches[nwa_column_check]],
                            axis=0,ignore_index=True)
    nwa_matches = nwa_matches.groupby('nwa_edge_id')['fraction_buffer'].sum().reset_index()
    nwa_roads = pd.merge(nwa_roads,nwa_matches,how='left',on=['nwa_edge_id']).fillna(0)
    nwa_roads[nwa_roads['fraction_buffer'] <= 0.5].to_file(store_intersections,layer='nwa_remaining_matches',driver='GPKG')

    nwa_edges = nwa_roads[nwa_roads['fraction_buffer'] <= 0.5].copy()
    match_edges = edges[~edges['id'].isin(all_matches)]
    road_classes = ['CLASS A','CLASS B','CLASS C','METRO']
    road_classes = ['OTHER','TRACK']
    match_edges = match_edges[~match_edges['_7'].isin(road_classes)]
    road_select = match_roads(nwa_edges,
                        match_edges,
                        geom_buffer=120,
                        fraction_intersection=0.9,
                        length_intersected=120,
                        save_buffer_file=store_intersections)
    road_select.to_file(store_intersections,layer='wider_buffer_match',driver='GPKG')

    # Get the road network id's and the NWA road ID's that uniquely match    
    df1, unique_matches = get_unique_matches(road_select)
    unique_matches.to_file(store_intersections,layer='unique_matches_wider_buffer',driver='GPKG')

    multiple_matches = road_select[road_select['id'].isin(df1[df1['count'] > 1]['id'].values.tolist())]
    multiple_matches['longest_section'] = multiple_matches['nwa_length']*multiple_matches['fraction_buffer']*multiple_matches['length_intersected']
    filter_multiple_matches = multiple_matches.groupby('id')['longest_section'].max().reset_index()
    filter_multiple_matches['longest_section_match'] = 1
    multiple_matches = pd.merge(multiple_matches,filter_multiple_matches,how='left',on=['id','longest_section']).fillna(0)
    multiple_matches[multiple_matches['longest_section_match'] == 1].to_file(store_intersections,
                                                layer='length_matches_filter_wider_buffer',driver='GPKG')


    network.edges = network.edges.rename(
        columns={
            "from_id": "from_node",
            "to_id": "to_node",
            "id": "edge_id",
            "_7": "road_class",
            "street_nam": "street_name",
            "street_typ": "street_type",
        }
    ).drop(
        columns=[
            "fnode_",
            "tnode_",
            "lpoly_",
            "rpoly_",
            "length",
            "shape_length",
            "road50west",
            "road50we_1",
        ]
    )
    # df.to_csv(os.path.join(processed_data_path,
    #         'networks','transport','roads_matches_count_2.csv'), index=False)

    # nwa_network = snkit.Network(edges=nwa_edges)
    # network = snkit.network.split_multilinestrings(network)
    # print ('* Done with splitting multilines')

    # network = snkit.network.snap_nodes(network)
    # print ('* Done with snapping nodes to edges')

    # network.nodes = snkit.network.drop_duplicate_geometries(network.nodes)
    # print ('* Done with dropping same geometries')
    # nwa_network = snkit.network.add_endpoints(nwa_network)
    # nwa_network = snkit.network.add_topology(nwa_network, id_col='id')

    # nwa_network.edges.to_file(os.path.join(processed_data_path,'networks','transport','nwa_roads.gpkg'), layer='edges', driver='GPKG')
    # nwa_network.nodes.to_file(os.path.join(processed_data_path,'networks','transport','nwa_roads.gpkg'), layer='nodes', driver='GPKG') 

    # nwa_edges = gpd.read_file(os.path.join(database_path,'nsdmb','GWP_Jamaica_NSP_Master_Geodatabase_v01.gdb'),
    #                     layer='roads_main_NWA')

    # create_network_from_nodes_and_edges(None,nwa_edges,
    #                                 'nwa_road',
    #                                 os.path.join(processed_data_path,'networks','transport',"nwa_roads_nomerge.gpkg"),by=None)

    # nwa_edges = gpd.read_file(os.path.join(processed_data_path,'networks','transport',"nwa_roads.gpkg"),layer='edges')
    # nwa_edges.set_crs(epsg=epsg_jamaica)

    # nwa_nodes = gpd.read_file(os.path.join(processed_data_path,'networks','transport',"nwa_roads.gpkg"),layer='nodes')
    # nwa_nodes.set_crs(epsg=epsg_jamaica) 

    # nwa_nodes['nid'] = nwa_nodes.geometry.progress_apply(
    #     lambda x: get_nearest_node(x, sindex_nodes,nodes, 'nid'))
    # nwa_nodes['geom'] = nwa_nodes.geometry.progress_apply(
    #     lambda x: get_nearest_node(x, sindex_nodes,nodes, 'geometry'))
    # nwa_nodes['node_distance'] = nwa_nodes.progress_apply(lambda x: x.geometry.distance(x.geom),axis=1)
    # nwa_nodes.drop(['geom'],axis=1,inplace=True)

    # nwa_nodes['eid'] = nwa_nodes.geometry.progress_apply(
    #     lambda x: get_nearest_node(x, sindex_edges,edges, 'eid'))
    # nwa_nodes['geom'] = nwa_nodes.geometry.progress_apply(
    #     lambda x: get_nearest_node(x, sindex_edges,edges, 'geometry'))
    # nwa_nodes['edge_distance'] = nwa_nodes.progress_apply(lambda x: x.geometry.distance(x.geom),axis=1)
    # nwa_nodes.drop(['geom'],axis=1,inplace=True)
    # nwa_nodes.to_csv(os.path.join(
    #                     processed_data_path,
    #                     'networks','transport',
    #                     'points_nearest.csv'),
    #                     index=False)

if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)