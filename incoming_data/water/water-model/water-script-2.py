import pandas as pd
import numpy as np
import geopandas as gpd
import math

# ########################################################################
# ########################################################################
# # potable
# ########################################################################
# ########################################################################

def roundup(x):
    return int(math.ceil(x / 50.0)) * 50

wsz_census = pd.read_csv('potable/processed/wsz_census.csv')
wsz_census = pd.read_csv('potable/processed/population_wsz.csv')
wsz_census['par'] = wsz_census['WSZONEID'].str[:3]

wsz_census_gb = wsz_census.groupby(['WSZONEID'])['population'].sum().to_frame().reset_index()
wsz_census_gb['par'] = wsz_census_gb['WSZONEID'].str[:3]
parish_names_df = pd.DataFrame({'par':['ann','eli','mar','cla','tho','ksa','tre','wes','cat','por','jam','man','han'],
                                'PARISH':['ST.ANN','ST.ELIZABETH','ST.MARY','CLARENDON','ST.THOMAS','KSA',
                                'TRELAWNY','WESTMORELAND','ST.CATHERINE','PORTLAND','ST.JAMES','MANCHESTER','HANOVER']})
wsz_census_gb = wsz_census_gb.merge(parish_names_df,on='par')

wsz_census['par'] = wsz_census['WSZONEID'].str[:3]
wsz_census_parish_pop = wsz_census.merge(parish_names_df,on='par')

parish_stats = pd.read_csv('potable/parish_stats.csv')
wsz_census_gb_par = wsz_census_gb.merge(parish_stats[['PARISH','% RESIDENTIAL CONNECTED 2010','Leakage rate 2010 ',
                                                      'Supply 2010','Population 2010','Population 2030','Commercial','Condominiums','Employees','Government','Residential','School']], 
                                                       on='PARISH', how='left').drop_duplicates()
wsz_census_gb_par['system_connected_pop_2010'] = wsz_census_gb_par['population']*wsz_census_gb_par['% RESIDENTIAL CONNECTED 2010']/100
wsz_census_gb_par['system_connected_pop_2010'] = wsz_census_gb_par['system_connected_pop_2010'] / 1000
wsz_census_gb_par['system_connected_pop_2010'] = wsz_census_gb_par['system_connected_pop_2010'].apply(np.ceil) * 1000

population_growth = pd.read_csv('potable/processed/pop_growth.csv')
population_growth_zone = population_growth.groupby(['WSZONEID'])['2020'].sum().to_frame().reset_index()
population_growth_zone['2011'] = population_growth.groupby(['WSZONEID'])['2011'].sum().to_frame().reset_index()['2011']
population_growth_zone['2030'] = population_growth.groupby(['WSZONEID'])['2030'].sum().to_frame().reset_index()['2030']
population_growth_zone['2040'] = population_growth.groupby(['WSZONEID'])['2040'].sum().to_frame().reset_index()['2040']
population_growth_zone['2050'] = population_growth.groupby(['WSZONEID'])['2050'].sum().to_frame().reset_index()['2050']
population_growth_zone['2011_2020_growth'] = (population_growth_zone['2020']-population_growth_zone['2011'])/population_growth_zone['2011']
population_growth_zone['2011_2030_growth'] = (population_growth_zone['2030']-population_growth_zone['2011'])/population_growth_zone['2011']
population_growth_zone['2011_2040_growth'] = (population_growth_zone['2040']-population_growth_zone['2011'])/population_growth_zone['2011']
population_growth_zone['2011_2050_growth'] = (population_growth_zone['2050']-population_growth_zone['2011'])/population_growth_zone['2011']

wsz_census_gb_par = pd.merge(wsz_census_gb_par, population_growth_zone, left_on='WSZONEID', right_on='WSZONEID')
wsz_census_gb_par['system_connected_pop_2020'] = wsz_census_gb_par['system_connected_pop_2010']+wsz_census_gb_par['system_connected_pop_2010']*wsz_census_gb_par['2011_2020_growth']
wsz_census_gb_par['system_connected_pop_2030'] = wsz_census_gb_par['system_connected_pop_2010']+wsz_census_gb_par['system_connected_pop_2010']*wsz_census_gb_par['2011_2030_growth']
wsz_census_gb_par['system_connected_pop_2040'] = wsz_census_gb_par['system_connected_pop_2010']+wsz_census_gb_par['system_connected_pop_2010']*wsz_census_gb_par['2011_2040_growth']
wsz_census_gb_par['system_connected_pop_2050'] = wsz_census_gb_par['system_connected_pop_2010']+wsz_census_gb_par['system_connected_pop_2010']*wsz_census_gb_par['2011_2050_growth']

GVA_per_sector_per_parish = pd.read_csv('potable/processed/GVA_per_zone_per_sector_per_parish_3.csv') # GVA in 1000 ## check GVA
GVA_per_sector_per_parish['building_type_nat_GVA'] = np.where(GVA_per_sector_per_parish['building_type'].isin(['Resort','Recreation']),'hotels/condominium',
                                                  np.where(GVA_per_sector_per_parish['building_type'].isin(['Commercial']),'commerical',
                                                  np.where(GVA_per_sector_per_parish['building_type'].isin(['Institutional']),'gov/school','nan')))
GVA_per_sector_per_parish = GVA_per_sector_per_parish[GVA_per_sector_per_parish['building_type']!='nan']

GVA_per_sector_per_parish['GVA (JD/day)'] = GVA_per_sector_per_parish['total_GDP']
GVA_per_sector_per_zone_gb = GVA_per_sector_per_parish.groupby(['WSZONEID'])['osm_id'].apply(list).to_frame().reset_index()
GVA_per_sector_per_zone_gb['GVA (JD/day)'] = GVA_per_sector_per_parish.groupby(['WSZONEID'])['GVA (JD/day)'].sum().to_frame().reset_index()['GVA (JD/day)']
GVA_per_sector_per_zone_gb['osm_id'] = GVA_per_sector_per_zone_gb['osm_id'].astype(str)

wsz_census_gb_par = pd.merge(wsz_census_gb_par, GVA_per_sector_per_zone_gb, left_on=['WSZONEID'], right_on=['WSZONEID'],how='left')
wsz_census_gb_par.to_csv('/soge-home/projects/mistral/jamaica-ccri/drought/calibrated_method/input/wsz_census_gb_par_update.csv')

Water_Supply_Zones = gpd.read_file('potable/Water_Supply_Zone.gpkg')
Water_Supply_Zones = Water_Supply_Zones.merge(wsz_census_gb_par[['WSZONEID','system_connected_pop_2010','system_connected_pop_2020','system_connected_pop_2030','system_connected_pop_2040','system_connected_pop_2050','GVA (JD/day)','osm_id']], on='WSZONEID') #
Water_Supply_Zones = gpd.GeoDataFrame(Water_Supply_Zones, crs="EPSG:3448",geometry=Water_Supply_Zones['geometry'])
Water_Supply_Zones.to_file('Water_Supply_Zone_update.gpkg', driver='GPKG')

parish_stats = pd.merge(parish_stats,wsz_census_gb[['PARISH','par']],on='PARISH', how='right').drop_duplicates()

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
potable_assets['asset_connected_GVA_other_sectors_JD/day'] = potable_assets['frac_pop']*potable_assets['GVA (JD/day)']
potable_assets['asset_connected_pop_2010'] = potable_assets['asset_connected_pop_2010']/1000
potable_assets['asset_connected_pop_2010'] = potable_assets['asset_connected_pop_2010'].apply(np.ceil) * 1000

potable_assets['par'] = potable_assets['WSZONEID'].str[:3]

# potable_assets = pd.merge(potable_assets, population_growth_zone, left_on='WSZONEID', right_on='WSZONEID')
potable_assets['asset_connected_pop_2020'] = potable_assets['asset_connected_pop_2010']+potable_assets['asset_connected_pop_2010']*potable_assets['2011_2020_growth']
potable_assets['asset_connected_pop_2030'] = potable_assets['asset_connected_pop_2010']+potable_assets['asset_connected_pop_2010']*potable_assets['2011_2030_growth']
potable_assets['asset_connected_pop_2040'] = potable_assets['asset_connected_pop_2010']+potable_assets['asset_connected_pop_2010']*potable_assets['2011_2040_growth']
potable_assets['asset_connected_pop_2050'] = potable_assets['asset_connected_pop_2010']+potable_assets['asset_connected_pop_2010']*potable_assets['2011_2050_growth']

potable_assets = pd.merge(potable_assets, parish_stats, left_on='par', right_on='par')

potable_facilities_NWC = gpd.read_file('/soge-home/projects/mistral/jamaica-ccri/processed_data/networks/water/ss/potable_facilities_NWC.gpkg')
potable_facilities_NWC = potable_facilities_NWC.drop(['asset_pop_new_2020'], axis=1)
potable_facilities_NWC_pop = potable_facilities_NWC.merge(potable_assets[['node_id','node_type','X','Y','WSZONEID','system_connected_pop_2010','asset_connected_pop_2010','asset_connected_pop_2020','asset_connected_pop_2030','asset_connected_pop_2040','asset_connected_pop_2050','asset_connected_GVA_other_sectors_JD/day','osm_id']], on='node_id', how='left').drop_duplicates()

potable_facilities_NWC_pop['W_GDP'] = 17648*(1e6/365)*0.8 * (potable_facilities_NWC_pop['asset_connected_pop_2010']/potable_facilities_NWC_pop['system_connected_pop_2010'])*(potable_facilities_NWC_pop['system_connected_pop_2010']/potable_assets.groupby('WSZONEID')['asset_connected_pop_2020'].max().to_frame().reset_index()['asset_connected_pop_2020'].sum())
print(17648*(1e6/365)*0.8,potable_facilities_NWC_pop.groupby('WSZONEID')['W_GDP'].max().to_frame().reset_index()['W_GDP'].sum(),potable_assets.groupby('WSZONEID')['asset_connected_pop_2020'].max().to_frame().reset_index()['asset_connected_pop_2020'].sum())

potable_facilities_NWC_economic = potable_facilities_NWC_pop[(potable_facilities_NWC_pop['node_type'].isin(['source','junction']))&(potable_facilities_NWC_pop['asset_connected_pop_2020']>0)]
potable_facilities_NWC_economic_m = potable_facilities_NWC_economic[['node_id','WSZONEID']].merge(GVA_per_sector_per_parish[['osm_id','A_GDP', 'B_GDP', 'C_GDP', 'D_GDP', 'E_GDP', 'F_GDP', 'G_GDP', 'H_GDP',
       'I_GDP', 'J_GDP', 'K_GDP', 'L_GDP', 'M_GDP', 'N_GDP', 'O_GDP',
       'total_GDP', 'GDP_unit','WSZONEID']], on='WSZONEID', how='left')
potable_facilities_NWC_economic_m = potable_facilities_NWC_economic_m[['node_id','osm_id','A_GDP', 'B_GDP', 'C_GDP', 'D_GDP', 'E_GDP', 'F_GDP', 'G_GDP', 'H_GDP',
       'I_GDP', 'J_GDP', 'K_GDP', 'L_GDP', 'M_GDP', 'N_GDP', 'O_GDP',
       'total_GDP', 'GDP_unit','WSZONEID']]
potable_facilities_NWC_economic_m.to_csv('potable_facilities_buildings_economic_activity_mapping.csv') # W_GDP
potable_facilities_NWC_economic_m = pd.read_csv('potable_facilities_buildings_economic_activity_mapping.csv')
print(potable_facilities_NWC_economic_m.columns)
potable_facilities_NWC_economic_m_2 = potable_facilities_NWC_pop[['node_id','W_GDP','WSZONEID']].merge(GVA_per_sector_per_parish[['osm_id','A_GDP', 'B_GDP', 'C_GDP', 'D_GDP', 'E_GDP', 'F_GDP', 'G_GDP', 'H_GDP',
       'I_GDP', 'J_GDP', 'K_GDP', 'L_GDP', 'M_GDP', 'N_GDP', 'O_GDP',
       'total_GDP', 'GDP_unit','WSZONEID']], on='WSZONEID', how='left')
potable_facilities_NWC_economic_m_2_gb = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['W_GDP'].mean().to_frame().reset_index()
potable_facilities_NWC_economic_m_2_gb['A_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['A_GDP'].sum().to_frame().reset_index()['A_GDP']
potable_facilities_NWC_economic_m_2_gb['B_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['B_GDP'].sum().to_frame().reset_index()['B_GDP']
potable_facilities_NWC_economic_m_2_gb['C_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['C_GDP'].sum().to_frame().reset_index()['C_GDP']
potable_facilities_NWC_economic_m_2_gb['D_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['D_GDP'].sum().to_frame().reset_index()['D_GDP']
potable_facilities_NWC_economic_m_2_gb['E_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['E_GDP'].sum().to_frame().reset_index()['E_GDP']
potable_facilities_NWC_economic_m_2_gb['F_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['F_GDP'].sum().to_frame().reset_index()['F_GDP']
potable_facilities_NWC_economic_m_2_gb['G_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['G_GDP'].sum().to_frame().reset_index()['G_GDP']
potable_facilities_NWC_economic_m_2_gb['H_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['H_GDP'].sum().to_frame().reset_index()['H_GDP']
potable_facilities_NWC_economic_m_2_gb['I_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['I_GDP'].sum().to_frame().reset_index()['I_GDP']
potable_facilities_NWC_economic_m_2_gb['J_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['J_GDP'].sum().to_frame().reset_index()['J_GDP']
potable_facilities_NWC_economic_m_2_gb['K_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['K_GDP'].sum().to_frame().reset_index()['K_GDP']
potable_facilities_NWC_economic_m_2_gb['L_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['L_GDP'].sum().to_frame().reset_index()['L_GDP']
potable_facilities_NWC_economic_m_2_gb['M_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['M_GDP'].sum().to_frame().reset_index()['M_GDP']
potable_facilities_NWC_economic_m_2_gb['N_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['N_GDP'].sum().to_frame().reset_index()['N_GDP']
potable_facilities_NWC_economic_m_2_gb['O_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['O_GDP'].sum().to_frame().reset_index()['O_GDP']
potable_facilities_NWC_economic_m_2_gb.to_csv('potable_facilities_dependent_economic_activity.csv')
# potable_facilities_NWC_economic_m_2_gb = pd.read_csv('potable_facilities_dependent_economic_activity.csv')
potable_facilities_NWC_pop['cost ($J) - lower bound'] = np.where(potable_facilities_NWC_pop['asset_type']=='well',5000000,potable_facilities_NWC_pop['cost ($J) - lower bound']) # same for irrigation
potable_facilities_NWC_pop['cost ($J) - upper bound'] = np.where(potable_facilities_NWC_pop['asset_type']=='well',20000000,potable_facilities_NWC_pop['cost ($J) - upper bound']) # same for irrigation
potable_facilities_NWC_pop['min_damage_cost'] = np.where(potable_facilities_NWC_pop['asset_type']=='well',5000000,potable_facilities_NWC_pop['min_damage_cost']) # same for irrigation
potable_facilities_NWC_pop['max_damage_cost'] = np.where(potable_facilities_NWC_pop['asset_type']=='well',20000000,potable_facilities_NWC_pop['max_damage_cost']) # same for irrigation

potable_facilities_NWC_pop.to_file('potable_facilities_NWC_pop_update.gpkg', driver='GPKG', layer='node')

pipelines_centroid = gpd.read_file('potable/processed/pipelines_centroid_filtered.gpkg')
pipelines_centroid['edge_id'] = [s.replace('pipeline', 'potable') for s in pipelines_centroid['edge_id'].to_list()] 
print(pipelines_centroid['edge_id'])
pipelines_NWC = gpd.read_file('/soge-home/projects/mistral/jamaica-ccri/processed_data/networks/water/pipelines_NWC_edges_pop_update.gpkg')
pipelines_NWC_pop = pipelines_NWC[['edge_id','geometry']].merge(pipelines_centroid[['edge_id','node_id']], on='edge_id', how='left').drop_duplicates()
print(pipelines_NWC_pop['node_id'])
pipelines_NWC_pop.to_file('pipelines_NWC_edges_pop_update.gpkg', driver='GPKG', layer='edge')
# pipelines_NWC_pop = gpd.read_file('pipelines_NWC_edges_pop_update.gpkg')
# print(pipelines_NWC_pop.columns,potable_facilities_NWC_economic_m.columns)
# print(pipelines_NWC_pop['node_id_x'],pipelines_NWC_pop['node_id_y'])
# pipelines_NWC_pop = pipelines_NWC_pop.rename(columns = {'node_id_x':'node_id'})
print(pipelines_NWC_pop['node_id'])
pipelines_NWC_pop = pipelines_NWC_pop.dropna(subset=['node_id'])
pipelines_NWC_pop = pipelines_NWC_pop.dropna(subset=['edge_id'])
potable_pipelines_NWC_economic_m = pd.merge(pipelines_NWC_pop[['edge_id','node_id']].head(50), potable_facilities_NWC_economic_m[['node_id','osm_id','A_GDP', 'B_GDP', 'C_GDP', 'D_GDP', 'E_GDP', 'F_GDP', 'G_GDP', 'H_GDP',
       'I_GDP', 'J_GDP', 'K_GDP', 'L_GDP', 'M_GDP', 'N_GDP', 'O_GDP',
       'total_GDP', 'GDP_unit','WSZONEID']], on='node_id', how='left')
       
print(potable_facilities_NWC_economic_m['node_id'])
print(potable_pipelines_NWC_economic_m['node_id'])
potable_pipelines_NWC_economic_m[['edge_id','osm_id','A_GDP', 'B_GDP', 'C_GDP', 'D_GDP', 'E_GDP', 'F_GDP', 'G_GDP', 'H_GDP',
       'I_GDP', 'J_GDP', 'K_GDP', 'L_GDP', 'M_GDP', 'N_GDP', 'O_GDP',
       'total_GDP', 'GDP_unit']].to_csv('/soge-home/projects/mistral/jamaica-ccri/processed_data/networks_economic_activity/potable_pipelines_buildings_economic_activity_mapping.csv') # W_GDP
potable_pipelines_NWC_economic_m = pd.merge(pipelines_NWC_pop[['edge_id','node_id']], potable_facilities_NWC_economic_m_2_gb[['node_id','W_GDP','A_GDP', 'B_GDP', 'C_GDP', 'D_GDP', 'E_GDP', 'F_GDP', 'G_GDP', 'H_GDP',
       'I_GDP', 'J_GDP', 'K_GDP', 'L_GDP', 'M_GDP', 'N_GDP', 'O_GDP']], on='node_id', how='left').to_csv('/soge-home/projects/mistral/jamaica-ccri/processed_data/networks_economic_activity/potable_pipelines_dependent_economic_activity.csv')
print(potable_pipelines_NWC_economic_m)
print(oliv)
########################################################################
########################################################################
# disruption analysis potable
########################################################################
########################################################################

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

NIC_IRRIGATION_SCHEMES = gpd.read_file('/soge-home/projects/mistral/jamaica-ccri/incoming_data/water/irrigation/raw/Irrigation Districts/NIC_IRRIGATION_SCHEMES.shp')

NIC_IRRIGATION_SCHEMES_convex = gpd.read_file('/soge-home/projects/mistral/jamaica-ccri/incoming_data/water/irrigation/raw/Irrigation Districts/NIC_IRRIGATION_SCHEMES_convex.shp')
NIC_IRRIGATION_SCHEMES_convex = NIC_IRRIGATION_SCHEMES_convex.rename(columns = {'DISTRICT':'District'})
NIC_IRRIGATION_SCHEMES_convex['District'] = np.where(NIC_IRRIGATION_SCHEMES_convex['District']=='BEACON','BEACON-LITTLE PARK',
                          np.where(NIC_IRRIGATION_SCHEMES_convex['District']=='MID - CLARENDON','MID CLARENDON',
                          np.where(NIC_IRRIGATION_SCHEMES_convex['District']=='RIO - COBRE','RIO COBRE', NIC_IRRIGATION_SCHEMES_convex['District'])))
NIC_IRRIGATION_SCHEMES_convex = NIC_IRRIGATION_SCHEMES_convex[~NIC_IRRIGATION_SCHEMES_convex.District.isin(NIC_IRRIGATION_SCHEMES['District'].drop_duplicates().to_list())]
NIC_IRRIGATION_SCHEMES = pd.concat([NIC_IRRIGATION_SCHEMES[['District','geometry']], NIC_IRRIGATION_SCHEMES_convex[['District','geometry']]])
NIC_IRRIGATION_SCHEMES['scheme_area_sqm'] = NIC_IRRIGATION_SCHEMES.geometry.area #* 1e-3 * 1e-3
NIC_IRRIGATION_SCHEMES['DXF_TEXT'] = np.where(NIC_IRRIGATION_SCHEMES['District']=='MID - CLARENDON', 'MID CLARENDON',
                          np.where(NIC_IRRIGATION_SCHEMES['District']=='RIO - COBRE', 'Rio Cobre', NIC_IRRIGATION_SCHEMES['District']))

agriculture = pd.read_csv('irrigation/processed/agriculture_gdp_intersect.csv')
# agriculture['intersect_area_sqm'] = agriculture['area'] #* 1e-3 * 1e-3
# agriculture['GDP_persqm'] = pd.to_numeric(agriculture['GDP_persqm'], errors='coerce')
# agriculture['GVA (JD/day)'] = agriculture['GDP_persqm']*agriculture['intersect_area_sqm']
agriculture['DXF_TEXT'] = np.where(agriculture['District']=='MID - CLARENDON', 'MID CLARENDON',
                          np.where(agriculture['District']=='RIO - COBRE', 'Rio Cobre', agriculture['District']))
agriculture_gb = agriculture.groupby('DXF_TEXT')['GVA (JD/day)'].sum().to_frame().reset_index()
agriculture_gb['GDP_persqm'] = agriculture.groupby('DXF_TEXT')['GDP_persqm'].mean().to_frame().reset_index()['GDP_persqm']

schemes = NIC_IRRIGATION_SCHEMES.merge(agriculture_gb, on='DXF_TEXT')
schemes_gb = schemes.groupby(['DXF_TEXT'])['GVA (JD/day)'].sum().to_frame().reset_index()
schemes_gb['GDP_persqm'] = schemes.groupby(['DXF_TEXT'])['GDP_persqm'].mean().to_frame().reset_index()['GDP_persqm']
schemes_gb['scheme_area_sqm'] = schemes.groupby(['DXF_TEXT'])['scheme_area_sqm'].mean().to_frame().reset_index()['scheme_area_sqm']
schemes_gb['GVA (JD/day)'] = schemes_gb['GDP_persqm']*schemes_gb['scheme_area_sqm']
NIC_IRRIGATION_SCHEMES = NIC_IRRIGATION_SCHEMES[['DXF_TEXT','geometry']].merge(schemes_gb, on='DXF_TEXT')
NIC_IRRIGATION_SCHEMES = gpd.GeoDataFrame(NIC_IRRIGATION_SCHEMES,crs="EPSG:3448",geometry=NIC_IRRIGATION_SCHEMES['geometry']) #.merge(schemes_gb, on='DXF_TEXT')
NIC_IRRIGATION_SCHEMES.drop_duplicates().to_file('NIC_IRRIGATION_SCHEMES_update_2.gpkg', driver='GPKG')

wells = pd.read_csv('irrigation/wells.csv')
wells['OBJECTID'] = np.linspace(0,wells.shape[0],wells.shape[0])
wells['node_id'] = 'well_'+wells['OBJECTID'].astype(float).astype(str)
wells_node_id = wells.to_csv('irrigation/wells_node_id.csv')

wells['DXF_TEXT'] = np.where(wells['DISTRICT']=='MID - CLARENDON', 'MID CLARENDON',
                    np.where(wells['DISTRICT']=='RIO - COBRE', 'Rio Cobre',
                    np.where(wells['DISTRICT']=='BEACON', 'BEACON-LITTLE PARK', wells['DISTRICT'])))
wells_gb = wells.groupby('DXF_TEXT')['node_id'].count().to_frame().reset_index()
wells_gb = wells_gb.rename(columns={'node_id':'no. wells'})

wells_merge = pd.merge(wells[['node_id','DXF_TEXT']],NIC_IRRIGATION_SCHEMES,on='DXF_TEXT',how='left').drop_duplicates()
wells_merge_2 = wells_merge.merge(wells_gb, on=['DXF_TEXT'],how='left').drop_duplicates()
wells_merge_2['sqm_per_node'] = wells_merge_2['scheme_area_sqm']/wells_merge_2['no. wells']
wells_merge_2['GDP/day_per_node'] = wells_merge_2['GVA (JD/day)']/wells_merge_2['no. wells']

irrigation_assets_NIC = gpd.read_file('/soge-home/projects/mistral/jamaica-ccri/processed_data/networks/water/ss/irrigation_assets_NIC.gpkg', layer ='nodes')
irrigation_assets_NIC['node_id'] = 'well_'+irrigation_assets_NIC['OBJECTID'].astype(float).astype(str)
irrigation_assets_NIC['cost ($J) - lower bound'] = np.where(irrigation_assets_NIC['asset_type']=='well',5000000,irrigation_assets_NIC['cost ($J) - lower bound']) # same for irrigation
irrigation_assets_NIC['cost ($J) - upper bound'] = np.where(irrigation_assets_NIC['asset_type']=='well',20000000,irrigation_assets_NIC['cost ($J) - upper bound']) # same for irrigation
irrigation_assets_NIC['min_damage_cost'] = np.where(irrigation_assets_NIC['asset_type']=='well',5000000,irrigation_assets_NIC['min_damage_cost']) # same for irrigation
irrigation_assets_NIC['max_damage_cost'] = np.where(irrigation_assets_NIC['asset_type']=='well',20000000,irrigation_assets_NIC['max_damage_cost']) # same for irrigation
irrigation_assets_NIC_are = irrigation_assets_NIC.merge(wells_merge_2[['node_id','DXF_TEXT', 'sqm_per_node','GDP/day_per_node']], on='node_id', how='left').drop_duplicates()
irrigation_assets_NIC_are.to_file('irrigation_assets_NIC_are_nodes_update.gpkg', driver='GPKG', layer='node')

# wells_NIC_economic_m = wells[['node_id','DXF_TEXT']].merge(agriculture[['fid','land_id','forest_id','GDP_persqm', 'sqm_per_node', 'GDP_unit', 'DXF_TEXT']], on='DXF_TEXT', how='left')
# wells_NIC_economic_m[['node_id','fid','land_id','forest_id','GDP_persqm', 'intersect_area_sqm', 'GDP_unit', 'GVA (JD/day)', 'DXF_TEXT']].to_csv('irrigation_nodes_agriculture_economic_activity_mapping.csv') # W_GDP

# wells_NIC_economic_m_2 = wells[['node_id','DXF_TEXT']].merge(agriculture[['fid','land_id','forest_id','GDP_persqm', 'intersect_area_sqm', 'GDP_unit', 'DXF_TEXT']], on='DXF_TEXT', how='left')
# wells_NIC_economic_m_2['GVA (JD/day)'] = wells_NIC_economic_m_2['GDP_persqm']*wells_NIC_economic_m_2['intersect_area_sqm']
# wells_NIC_economic_m_3 = wells_NIC_economic_m_2.groupby('node_id')['GVA (JD/day)'].sum().to_frame().reset_index()  
# wells_NIC_economic_m_3.to_csv('irrigation_nodes_dependent_economic_activity.csv')

pipelines_centroids = gpd.read_file('irrigation/processed/pipe_centroids_output_from_qgis_2.gpkg') ## read from qgis
pipelines_centroids['edge_id'] = 'pipeline_'+pipelines_centroids['OBJECTID'].astype(float).astype(str) 
pipelines_centroids = pipelines_centroids.merge(irrigation_assets_NIC_are[['DXF_TEXT','node_id','sqm_per_node','GDP/day_per_node']], on='node_id')
    
canals_centroids = gpd.read_file('irrigation/processed/canal_centroids_output_from_qgis_2.gpkg') ## read from qgis
canals_centroids['edge_id'] = 'pipeline_'+canals_centroids['OBJECTID'].astype(float).astype(str) 
canals_centroids = canals_centroids.merge(irrigation_assets_NIC_are[['DXF_TEXT','node_id','sqm_per_node','GDP/day_per_node']], on='node_id')

edges_centroids = pd.concat([pipelines_centroids,canals_centroids])

irrigation_assets_NIC_edges = gpd.read_file('irrigation/irrigation_edges.gpkg')
irrigation_assets_NIC_edges['edge_id'] = 'pipeline_'+irrigation_assets_NIC_edges['OBJECTID'].astype(str)
irrigation_assets_NIC_edges['edge_id'] = irrigation_assets_NIC_edges['edge_id'].astype(str)

edges_centroids['edge_id'] = edges_centroids['edge_id'].astype(str)

irrigation_assets_NIC_are_egdes = irrigation_assets_NIC_edges.merge(edges_centroids[['edge_id','node_id','sqm_per_node','GDP/day_per_node']], on='edge_id', how='left').drop_duplicates()
irrigation_assets_NIC_are_egdes = gpd.GeoDataFrame(irrigation_assets_NIC_are_egdes,crs="EPSG:3448",geometry=irrigation_assets_NIC_are_egdes['geometry'])
irrigation_assets_NIC_are_egdes.to_file('irrigation_assets_NIC_are_edges_update.gpkg', driver='GPKG', layer='edge') #.drop_duplicates()

# irrigation_edges_NIC_economic_m = edges_centroids[['edge_id','node_id']].merge(irrigation_assets_NIC_are[['node_id','fid','land_id','forest_id', 'GVA (JD/day)', 'DXF_TEXT']], on='node_id', how='left')
# irrigation_edges_NIC_economic_m.to_csv('irrigation_edges_agriculture_economic_activity_mapping.csv') # W_GDP
# irrigation_edges_NIC_economic_m = edges_centroids[['edge_id','node_id']].merge(irrigation_assets_NIC_are, on='node_id', how='left')
# irrigation_edges_NIC_economic_m.to_csv('irrigation_edges_dependent_economic_activity.csv')
print(oliv)

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
    
