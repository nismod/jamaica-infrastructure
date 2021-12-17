import pandas as pd
import numpy as np
import geopandas as gpd
import math

#########################################################################
#########################################################################
# potable
#########################################################################
#########################################################################

def roundup(x):
    return int(math.ceil(x / 50.0)) * 50

wsz_census = pd.read_csv('potable/processed/wsz_census.csv')
# wsz_census['PARISH'] = np.where(wsz_census['PARISH']=='ST.ANDREW', 'KSA',
                          # np.where(wsz_census['PARISH']=='KINGSTON', 'KSA', wsz_census['PARISH']) )
wsz_census_gb = wsz_census.groupby(['WSZONEID'])['TOTAL_POP'].sum().to_frame().reset_index()
wsz_census_gb['par'] = wsz_census_gb['WSZONEID'].str[:3]

parish_names_df = pd.DataFrame({'par':['ann','eli','mar','cla','tho','ksa','tre','wes','cat','por','jam'],
                                'PARISH':['ST.ANN','ST.ELIZABETH','ST.MARY','CLARENDON','ST.THOMAS','KSA',
                                'TRELAWNY','WESTMORELAND','ST.CATHERINE','PORTLAND','ST.JAMES']})
wsz_census_gb = wsz_census_gb.merge(parish_names_df,on='par')

parish_stats = pd.read_csv('potable/parish_stats.csv')
wsz_census_gb_par = wsz_census_gb.merge(parish_stats[['PARISH','% RESIDENTIAL CONNECTED 2010','Leakage rate 2010 ',
                                                      'Supply 2010','Population 2010','Population 2030']], 
                                                       on='PARISH', how='left').drop_duplicates()
wsz_census_gb_par['Population 2030'] = pd.to_numeric(wsz_census_gb_par['Population 2030'], errors='coerce')
wsz_census_gb_par['Population 2010'] = pd.to_numeric(wsz_census_gb_par['Population 2010'], errors='coerce')
wsz_census_gb_par['Growth'] = wsz_census_gb_par['Population 2030']/wsz_census_gb_par['Population 2010']
wsz_census_gb_par['system_connected_pop_2010'] = wsz_census_gb_par['TOTAL_POP']*wsz_census_gb_par['% RESIDENTIAL CONNECTED 2010']/100
wsz_census_gb_par['system_connected_pop_2030'] = wsz_census_gb_par['system_connected_pop_2010']*wsz_census_gb_par['Growth']

# check = wsz_census_gb_par.groupby(['WSZONEID'])['system_connected_pop_2010'].count().to_frame().reset_index()
# print(check[check['system_connected_pop_2010']==1])

wsz_potable_assets = pd.read_csv('potable/processed/wsz_potable_supply_assets_3.csv')
wsz_potable_assets['OBJECTID'] = pd.to_numeric(wsz_potable_assets['OBJECTID'])
wsz_potable_assets['node_id'] = wsz_potable_assets['Type'].astype(str)+'_'+wsz_potable_assets['OBJECTID'].astype(float).astype(str)

wsz_potable_assets_merge = wsz_potable_assets.merge(wsz_census_gb_par, on='WSZONEID',how='left').drop_duplicates()
wsz_potable_assets_merge['Type'] = np.where(wsz_potable_assets_merge['LOCATION'].isin(['Mona Dam', 'Hermitage Dam']), 'dam', wsz_potable_assets_merge['Type'])

source = ['Spring','Treatment Plant','Intake','Production Well','River Source','Filter Plant','Entombment', 'dam']
junction = ['Relift Station','Pump Station','Booster Station']
end = ['Storage Tank','Reservoir']
wsz_potable_assets_merge['node_type'] = np.where(wsz_potable_assets_merge['Type'].isin(source),'source',
                                        np.where(wsz_potable_assets_merge['Type'].isin(junction),'junction',
                                        np.where(wsz_potable_assets_merge['Type'].isin(end),'end','nan')))
wsz_potable_assets_merge_gb = wsz_potable_assets_merge.groupby(['WSZONEID','node_type'])['LOCATION'].count().to_frame().reset_index()
wsz_potable_assets_merge_gb_end = wsz_potable_assets_merge_gb[wsz_potable_assets_merge_gb['node_type']=='end']
wsz_potable_assets_merge_gb_end = wsz_potable_assets_merge_gb_end.rename(columns = {'LOCATION':'no. ends'})

potable_assets = wsz_potable_assets_merge.merge(wsz_potable_assets_merge_gb_end[['WSZONEID','no. ends']], on='WSZONEID',how='left').drop_duplicates()
potable_assets['frac_pop'] = np.where(potable_assets['Type'].isin(source),1,
                                        np.where(potable_assets['Type'].isin(junction),1,
                                        np.where(potable_assets['Type'].isin(end),1/potable_assets['no. ends'],0)))
potable_assets['asset_connected_pop_2010'] = potable_assets['frac_pop']*potable_assets['system_connected_pop_2010']
potable_assets['asset_connected_pop_2030'] = potable_assets['frac_pop']*potable_assets['system_connected_pop_2030']

# check = potable_assets.groupby(['node_id','WSZONEID'])['asset_connected_pop_2010'].count().to_frame().reset_index()
# print(check[check['asset_connected_pop_2010']>1])
# print(oliv)

potable_facilities_NWC = gpd.read_file('/soge-home/projects/mistral/jamaica-ccri/processed_data/networks/water/potable_facilities_NWC.gpkg')
potable_facilities_NWC = potable_facilities_NWC.drop(['asset_pop_new_2010','asset_pop_new_2020','asset_pop_new_2030'], axis=1)
potable_facilities_NWC_pop = potable_facilities_NWC.merge(potable_assets[['node_id','WSZONEID','asset_connected_pop_2010']], on='node_id', how='left').drop_duplicates()
#potable_facilities_NWC_pop['asset_connected_pop_2010'] = potable_facilities_NWC_pop['asset_connected_pop_2010'].apply(roundup)
potable_facilities_NWC_pop.to_file('potable_facilities_NWC_pop.gpkg', driver='GPKG')

pipelines_centroid = pd.read_csv('potable/processed/pipelines_centroid.csv')
pipelines_centroid['OBJECTID'] = pd.to_numeric(pipelines_centroid['OBJECTID'])
pipelines_centroid['edge_id'] = 'potable_'+pipelines_centroid['OBJECTID'].astype(float).astype(str)
for it,row in pipelines_centroid.iterrows(): #canal/pipe centroid
    pipeline_id = row.OBJECTID #OBJECTID
    pipelines_Y = row.Y
    pipelines_X = row.X
    potable_assets['lat_diff'] = abs(potable_assets['Y'] - pipelines_Y)
    potable_assets['lon_diff'] = abs(potable_assets['X'] - pipelines_X)
    potable_assets['sqrt'] = np.sqrt(np.square(potable_assets['lat_diff'])+np.square(potable_assets['lat_diff']))
    closest_wrz = potable_assets[potable_assets['sqrt']==potable_assets['sqrt'].min()]['WSZONEID'].values[0]
    asset_connected_pop_2010 = potable_assets[potable_assets['sqrt']==potable_assets['sqrt'].min()]['asset_connected_pop_2010'].values[0]
    asset_connected_pop_2030 = potable_assets[potable_assets['sqrt']==potable_assets['sqrt'].min()]['asset_connected_pop_2030'].values[0]
    pipelines_centroid.loc[it,'WSZONEID'] = closest_wrz
    pipelines_centroid.loc[it,'asset_connected_pop_2010'] = asset_connected_pop_2010
    pipelines_centroid.loc[it,'asset_connected_pop_2030'] = asset_connected_pop_2030

pipelines_NWC = gpd.read_file('/soge-home/projects/mistral/jamaica-ccri/processed_data/networks/water/pipelines_NWC.gpkg')
pipelines_NWC_pop = pipelines_NWC.merge(pipelines_centroid[['edge_id','WSZONEID','asset_connected_pop_2010','asset_connected_pop_2030']], on='edge_id', how='left').drop_duplicates()
#pipelines_NWC_pop = pipelines_NWC_pop.dropna(subset=['asset_connected_pop_2010'])
#pipelines_NWC_pop['asset_connected_pop_2010'] = pipelines_NWC_pop['asset_connected_pop_2010'].apply(roundup)
pipelines_NWC_pop.to_file('pipelines_NWC_edges_pop.gpkg', driver='GPKG')

#########################################################################
#########################################################################
# disruption analysis potable
#########################################################################
#########################################################################

# damaged_nodes_df = pd.DataFrame([])
# merged_nodes = pd.merge(potable_assets,damaged_nodes_df,on='node_id')
# pop_disrupted_nodes = merged_nodes.groupby('WSZONEID')['asset_connected_pop_2010'].max().reset_index().to_frame() ### sink assets are storage tanks that are assumed not to be vulnerable

# damaged_edges_df = pd.DataFrame([])
# merged_edges = pd.merge(pipelines_centroid,damaged_edges_df,on='node_id')
# merged_edges_sink = merged_edges[merged_edges['closest_node_type']=='sink'] 
# merged_edges_sink_gb = merged_edges_sink.groupby(['closest_node_id','WSZONEID'])['asset_connected_pop_2010'].max().reset_index().to_frame()
# merged_edges_sink_gb_2 = merged_edges_sink_gb.groupby('WSZONEID')['asset_connected_pop_2010'].sum().reset_index().to_frame()
# merged_edges_concat = pd.concat([merged_edges_sink_gb_2, merged_edges_gb])
# pop_disrupted_ = merged_edges_concat.groupby('WSZONEID')['asset_connected_pop_2010'].max().reset_index().to_frame()

#########################################################################
#########################################################################
# irrigation
#########################################################################
#########################################################################

schemes = pd.read_csv('irrigation/schemes.csv')
wells = pd.read_csv('irrigation/wells.csv')
wells['OBJECTID'] = np.linspace(0,wells.shape[0],wells.shape[0])
wells['node_id'] = 'well_'+wells['OBJECTID'].astype(float).astype(str)

scheme_centroids = pd.read_csv('irrigation/processed/schemes_centroids.csv')
pipelines_centroids = pd.read_csv('irrigation/processed/pipe_centroids.csv')
pipelines_centroids['edge_id'] = 'pipeline_'+pipelines_centroids['OBJECTID'].astype(float).astype(str) 

canals_centroids = pd.read_csv('irrigation/processed/canal_centroids.csv')

wells['DISTRICT'].drop_duplicates()
wells['DXF_TEXT'] = np.where(wells['DISTRICT']=='MID - CLARENDON', 'MID CLARENDON',
                    np.where(wells['DISTRICT']=='RIO - COBRE', 'Rio Cobre', wells['DISTRICT']))

wells_merge = pd.merge(wells, schemes,on='DXF_TEXT',how='left').drop_duplicates()
wells_merge_gb = wells_merge.groupby(['DXF_TEXT','ACRES'])['NAME'].count().to_frame().reset_index()
wells_merge_gb = wells_merge_gb.rename(columns={'NAME':'no. wells'})
wells_merge_2 = wells_merge.merge(wells_merge_gb, on=['DXF_TEXT','ACRES'],how='left').drop_duplicates()
wells_merge_2['acres_per_node'] = wells_merge_2['ACRES']/wells_merge_2['no. wells']
#wells_merge_2['acres_per_node'] = wells_merge_2['acres_per_node'].apply(roundup)
 
irrigation_assets_NIC = gpd.read_file('/soge-home/projects/mistral/jamaica-ccri/processed_data/networks/water/irrigation_assets_NIC.gpkg', layer ='nodes')
irrigation_assets_NIC_are = irrigation_assets_NIC.merge(wells_merge_2[['node_id','acres_per_node']], on='node_id', how='left').drop_duplicates()
irrigation_assets_NIC_are.to_file('irrigation_assets_NIC_are_nodes.gpkg', driver='GPKG')

for it,row in pipelines_centroids.iterrows():
    pipeline_id = row.OBJECTID 
    pipelines_Y = row.Y
    pipelines_X = row.X
    pipelines_area = row.Size
    scheme_centroids['lat_diff'] = abs(scheme_centroids['Y'] - pipelines_Y)
    scheme_centroids['lon_diff'] = abs(scheme_centroids['X'] - pipelines_X)
    scheme_centroids['sqrt'] = np.sqrt(np.square(scheme_centroids['lat_diff'])+np.square(scheme_centroids['lat_diff']))
    closest_scheme = scheme_centroids[scheme_centroids['sqrt']==scheme_centroids['sqrt'].min()]['DXF_TEXT'].values[0]
    pipelines_centroids.loc[it,'DXF_TEXT'] = closest_scheme

pipelines_centroids_merge = pd.merge(pipelines_centroids,schemes[['DXF_TEXT','ACRES']],on='DXF_TEXT',how='left').drop_duplicates()
pipelines_centroids_merge_gb = pipelines_centroids_merge.groupby(['DXF_TEXT','ACRES'])['Size'].sum().to_frame().reset_index()
pipelines_centroids_merge_gb = pipelines_centroids_merge_gb.rename(columns={'Size':'sum of pipe size'})
pipelines_centroids_merge_2 = pd.merge(pipelines_centroids_merge, pipelines_centroids_merge_gb, on=['DXF_TEXT','ACRES'], how='left').drop_duplicates()
pipelines_centroids_merge_2['acres_per_edge'] = pipelines_centroids_merge_2['ACRES']*pipelines_centroids_merge_2['Size']/pipelines_centroids_merge_2['sum of pipe size']

for it,row in canals_centroids.iterrows(): #canal/pipe centroid
    pipeline_id = row.OBJECTID #
    pipelines_Y = row.Y
    pipelines_X = row.X
    scheme_centroids['lat_diff'] = abs(scheme_centroids['Y'] - pipelines_Y)
    scheme_centroids['lon_diff'] = abs(scheme_centroids['X'] - pipelines_X)
    scheme_centroids['sqrt'] = np.sqrt(np.square(scheme_centroids['lat_diff'])+np.square(scheme_centroids['lat_diff']))
    closest_scheme = scheme_centroids[scheme_centroids['sqrt']==scheme_centroids['sqrt'].min()]['DXF_TEXT'].values[0]
    canals_centroids.loc[it,'DXF_TEXT'] = closest_scheme
canals_centroids['DXF_TEXT'] = np.where(canals_centroids['Parish']=='Clarendon','MID CLARENDON',canals_centroids['DXF_TEXT'])

canals_centroids_merge = pd.merge(canals_centroids,schemes,on='DXF_TEXT',how='left').drop_duplicates()
canals_centroids_merge_gb = canals_centroids_merge.groupby(['DXF_TEXT','ACRES'])['LENGTH'].sum().to_frame().reset_index()
canals_centroids_merge_gb = canals_centroids_merge_gb.rename(columns={'LENGTH':'sum of canal length'})
canals_centroids_merge_2 = canals_centroids_merge.merge(canals_centroids_merge_gb, on=['DXF_TEXT','ACRES'],how='left').drop_duplicates()
canals_centroids_merge_2['acres_per_edge'] = canals_centroids_merge_2['ACRES']*canals_centroids_merge_2['LENGTH']/canals_centroids_merge_2['sum of canal length']

canals_centroids_merge_2['acres_per_edge'] = np.where(canals_centroids_merge_2['DXF_TEXT']=='MID CLARENDON',canals_centroids_merge_2['AREA_SERVE'],canals_centroids_merge_2['acres_per_edge'])
edges_centroids = pd.concat([pipelines_centroids_merge_2,canals_centroids_merge_2])
#edges_centroids = edges_centroids.dropna(subset=['acres_per_edge'])
#edges_centroids['acres_per_edge'] = edges_centroids['acres_per_edge'].apply(roundup)

irrigation_assets_NIC_edges = gpd.read_file('/soge-home/projects/mistral/jamaica-ccri/processed_data/networks/water/irrigation_assets_NIC.gpkg', layer ='edges')
irrigation_assets_NIC_are_egdes = irrigation_assets_NIC.merge(edges_centroids[['edge_id','acres_per_edge','DXF_TEXT']], on='edge_id', how='left').drop_duplicates()

# irrigation_assets_NIC_are_egdes.to_file('irrigation_assets_NIC_are_edges.gpkg', driver='GPKG').drop_duplicates()

#########################################################################
#########################################################################
# disruption analysis irrigation
#########################################################################
#########################################################################

# damaged_nodes_df = pd.DataFrame([])
# merged_nodes = pd.merge(wells_merge_2,damaged_nodes_df,on='node_id')
# merged_nodes_gb = merged_nodes.groupby('DXF_TEXT')['acres_per_node'].sum().reset_index().to_frame()

# damaged_edges_df = pd.DataFrame([])
# merged_edges = pd.merge(edges_centroids,damaged_edges_df,on='edge_id')
# merged_edges_gb = merged_edges.groupby('DXF_TEXT')['acres_per_edge'].sum().reset_index().to_frame()
    
