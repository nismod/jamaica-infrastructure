import pandas as pd
import numpy as np
import geopandas as gpd
import math


########################################################################
######################################################################## checks
########################################################################

# check system connected == asset
# wsz_census_gb_par = pd.read_csv('/soge-home/projects/mistral/jamaica-ccri/drought/calibrated_method/input/wsz_census_gb_par_update.csv')
# print(wsz_census_gb_par.columns)
# print(wsz_census_gb_par['system_connected_pop_2010'].sum(),wsz_census_gb_par['system_connected_pop_2020'].sum())

# potable_facilities_NWC_pop = gpd.read_file('potable_facilities_NWC_pop_update.gpkg')
# print(potable_facilities_NWC_pop.groupby('WSZONEID')['asset_connected_pop_2010'].max().to_frame().reset_index()['asset_connected_pop_2010'].sum())

# GVA_per_sector_per_parish = pd.read_csv('potable/processed/GVA_per_zone_per_sector_per_parish_3.csv') # GVA in 1000 ## check GVA
# GVA_per_sector_per_parish['building_type_nat_GVA'] = np.where(GVA_per_sector_per_parish['building_type'].isin(['Resort','Recreation']),'hotels/condominium',
                                                  # np.where(GVA_per_sector_per_parish['building_type'].isin(['Commercial']),'commerical',
                                                  # np.where(GVA_per_sector_per_parish['building_type'].isin(['Institutional']),'gov/school','nan')))
# GVA_per_sector_per_parish = GVA_per_sector_per_parish[GVA_per_sector_per_parish['building_type']!='nan']
# print('A_GDP', GVA_per_sector_per_parish['A_GDP'].sum())
# print('osm_id',GVA_per_sector_per_parish['osm_id'].count())

# potable_facilities_NWC_economic_m = pd.read_csv('potable_facilities_buildings_economic_activity_mapping.csv') # W_GDP
# potable_facilities_NWC_economic_m_wsz = potable_facilities_NWC_economic_m.groupby(['WSZONEID','node_id'])['osm_id'].count().to_frame().reset_index()#['osm_id'].sum())
# potable_facilities_NWC_economic_m_wsz['A_GDP'] = potable_facilities_NWC_economic_m.groupby(['WSZONEID','node_id'])['A_GDP'].sum().to_frame().reset_index()['A_GDP']
# print('osm_id',potable_facilities_NWC_economic_m_wsz.groupby(['WSZONEID'])['osm_id'].mean().to_frame().reset_index()['osm_id'].sum())
# print('A_GDP',potable_facilities_NWC_economic_m_wsz.groupby(['WSZONEID'])['A_GDP'].mean().to_frame().reset_index()['A_GDP'].sum())

# potable_facilities_NWC_economic_m_2_gb = pd.read_csv('potable_facilities_dependent_economic_activity.csv')
# print('W_GDP',potable_facilities_NWC_economic_m_2_gb.groupby('WSZONEID')['W_GDP'].max().to_frame().reset_index()['W_GDP'].sum(), 17648*(1e6/365)*0.8)
# print('A_GDP',potable_facilities_NWC_economic_m_2_gb.groupby('WSZONEID')['A_GDP'].max().to_frame().reset_index()['A_GDP'].sum())

# # --> same for are and gdp in agriculture areas
# # check system area == well area sum
# NIC_IRRIGATION_SCHEMES = gpd.read_file('NIC_IRRIGATION_SCHEMES_update.gpkg')
# print(NIC_IRRIGATION_SCHEMES.columns)
# print('area',NIC_IRRIGATION_SCHEMES['scheme_area_sqm'].sum())
# print('area per scheme', NIC_IRRIGATION_SCHEMES[['DXF_TEXT','scheme_area_sqm']])

# irrigation_assets_NIC_are = gpd.read_file('irrigation_assets_NIC_are_nodes_update.gpkg')
# print(irrigation_assets_NIC_are.columns)
# print('area',irrigation_assets_NIC_are['sqm_per_node'])
# print('gdp',irrigation_assets_NIC_are['GDP/day_per_node'].sum())

# agriculture = pd.read_csv('irrigation/processed/agriculture_gdp_intersect.csv')
# agriculture['intersect_area_sqm'] = agriculture['area'] #* 1e-3 * 1e-3
# agriculture['GDP_persqm'] = pd.to_numeric(agriculture['GDP_persqm'], errors='coerce')
# agriculture['GVA (JD/day)'] = agriculture['GDP_persqm']*agriculture['intersect_area_sqm']
# agriculture['DXF_TEXT'] = np.where(agriculture['District']=='MID - CLARENDON', 'MID CLARENDON',
                          # np.where(agriculture['District']=='RIO - COBRE', 'Rio Cobre', agriculture['District']))
# agriculture_gb = agriculture.groupby('DXF_TEXT')['GVA (JD/day)'].sum().to_frame().reset_index()
# agriculture_gb['GDP_persqm'] = agriculture.groupby('DXF_TEXT')['GDP_persqm'].mean().to_frame().reset_index()['GDP_persqm']
# agriculture_gb['economic_interection_area'] = agriculture.groupby('DXF_TEXT')['intersect_area_sqm'].sum().to_frame().reset_index()['intersect_area_sqm']
# print('area per scheme', agriculture_gb[['DXF_TEXT','economic_interection_area','GVA (JD/day)']])

# print('area',agriculture_gb['economic_interection_area'].sum())
# print('gdp',agriculture['GVA (JD/day)'].sum())

# print(oliv)

########################################################################
########################################################################
########################################################################

# ########################################################################
# ########################################################################
# # potable
# ########################################################################
# ########################################################################

# def roundup(x):
    # return int(math.ceil(x / 50.0)) * 50

# wsz_census = pd.read_csv('potable/processed/wsz_census.csv')
# wsz_census = pd.read_csv('potable/processed/population_wsz.csv')
# wsz_census['par'] = wsz_census['WSZONEID'].str[:3]

# wsz_census_gb = wsz_census.groupby(['WSZONEID'])['population'].sum().to_frame().reset_index()
# wsz_census_gb['par'] = wsz_census_gb['WSZONEID'].str[:3]
# parish_names_df = pd.DataFrame({'par':['ann','eli','mar','cla','tho','ksa','tre','wes','cat','por','jam','man','han'],
                                # 'PARISH':['ST.ANN','ST.ELIZABETH','ST.MARY','CLARENDON','ST.THOMAS','KSA',
                                # 'TRELAWNY','WESTMORELAND','ST.CATHERINE','PORTLAND','ST.JAMES','MANCHESTER','HANOVER']})
# wsz_census_gb = wsz_census_gb.merge(parish_names_df,on='par')

# wsz_census['par'] = wsz_census['WSZONEID'].str[:3]
# wsz_census_parish_pop = wsz_census.merge(parish_names_df,on='par')

# parish_stats = pd.read_csv('potable/parish_stats.csv')
# wsz_census_gb_par = wsz_census_gb.merge(parish_stats[['PARISH','% RESIDENTIAL CONNECTED 2010','Leakage rate 2010 ',
                                                      # 'Supply 2010','Population 2010','Population 2030','Commercial','Condominiums','Employees','Government','Residential','School']], 
                                                       # on='PARISH', how='left').drop_duplicates()
# # # # # wsz_census_gb_par['Population 2030'] = pd.to_numeric(wsz_census_gb_par['Population 2030'], errors='coerce')
# # # # # wsz_census_gb_par['Population 2010'] = pd.to_numeric(wsz_census_gb_par['Population 2010'], errors='coerce')
# # # # # wsz_census_gb_par['Growth'] = (wsz_census_gb_par['Population 2030']-wsz_census_gb_par['Population 2010'])/wsz_census_gb_par['Population 2010']
# wsz_census_gb_par['system_connected_pop_2010'] = wsz_census_gb_par['population']*wsz_census_gb_par['% RESIDENTIAL CONNECTED 2010']/100
# # # # # wsz_census_gb_par['system_connected_pop_2010'] = (wsz_census_gb_par['system_connected_pop_2010']+wsz_census_gb_par['system_connected_pop_2010']*wsz_census_gb_par['Growth'])/1000
# wsz_census_gb_par['system_connected_pop_2010'] = wsz_census_gb_par['system_connected_pop_2010'] / 1000
# wsz_census_gb_par['system_connected_pop_2010'] = wsz_census_gb_par['system_connected_pop_2010'].apply(np.ceil) * 1000

# population_growth = pd.read_csv('potable/processed/pop_growth.csv')
# population_growth_zone = population_growth.groupby(['WSZONEID'])['2020'].sum().to_frame().reset_index()
# population_growth_zone['2011'] = population_growth.groupby(['WSZONEID'])['2011'].sum().to_frame().reset_index()['2011']
# population_growth_zone['2030'] = population_growth.groupby(['WSZONEID'])['2030'].sum().to_frame().reset_index()['2030']
# population_growth_zone['2040'] = population_growth.groupby(['WSZONEID'])['2040'].sum().to_frame().reset_index()['2040']
# population_growth_zone['2050'] = population_growth.groupby(['WSZONEID'])['2050'].sum().to_frame().reset_index()['2050']
# population_growth_zone['2011_2020_growth'] = (population_growth_zone['2020']-population_growth_zone['2011'])/population_growth_zone['2011']
# population_growth_zone['2011_2030_growth'] = (population_growth_zone['2030']-population_growth_zone['2011'])/population_growth_zone['2011']
# population_growth_zone['2011_2040_growth'] = (population_growth_zone['2040']-population_growth_zone['2011'])/population_growth_zone['2011']
# population_growth_zone['2011_2050_growth'] = (population_growth_zone['2050']-population_growth_zone['2011'])/population_growth_zone['2011']

# wsz_census_gb_par = pd.merge(wsz_census_gb_par, population_growth_zone, left_on='WSZONEID', right_on='WSZONEID')
# wsz_census_gb_par['system_connected_pop_2020'] = wsz_census_gb_par['system_connected_pop_2010']+wsz_census_gb_par['system_connected_pop_2010']*wsz_census_gb_par['2011_2020_growth']
# wsz_census_gb_par['system_connected_pop_2030'] = wsz_census_gb_par['system_connected_pop_2010']+wsz_census_gb_par['system_connected_pop_2010']*wsz_census_gb_par['2011_2030_growth']
# wsz_census_gb_par['system_connected_pop_2040'] = wsz_census_gb_par['system_connected_pop_2010']+wsz_census_gb_par['system_connected_pop_2010']*wsz_census_gb_par['2011_2040_growth']
# wsz_census_gb_par['system_connected_pop_2050'] = wsz_census_gb_par['system_connected_pop_2010']+wsz_census_gb_par['system_connected_pop_2010']*wsz_census_gb_par['2011_2050_growth']

# GVA_per_sector_per_parish = pd.read_csv('potable/processed/GVA_per_zone_per_sector_per_parish_3.csv') # GVA in 1000 ## check GVA
# GVA_per_sector_per_parish['building_type_nat_GVA'] = np.where(GVA_per_sector_per_parish['building_type'].isin(['Resort','Recreation']),'hotels/condominium',
                                                  # np.where(GVA_per_sector_per_parish['building_type'].isin(['Commercial']),'commerical',
                                                  # np.where(GVA_per_sector_per_parish['building_type'].isin(['Institutional']),'gov/school','nan')))
# GVA_per_sector_per_parish = GVA_per_sector_per_parish[GVA_per_sector_per_parish['building_type']!='nan']

# GVA_per_sector_per_parish['GVA (JD/day)'] = GVA_per_sector_per_parish['total_GDP']
# GVA_per_sector_per_zone_gb = GVA_per_sector_per_parish.groupby(['WSZONEID'])['osm_id'].apply(list).to_frame().reset_index()
# GVA_per_sector_per_zone_gb['GVA (JD/day)'] = GVA_per_sector_per_parish.groupby(['WSZONEID'])['GVA (JD/day)'].sum().to_frame().reset_index()['GVA (JD/day)']
# GVA_per_sector_per_zone_gb['osm_id'] = GVA_per_sector_per_zone_gb['osm_id'].astype(str)

# wsz_census_gb_par = pd.merge(wsz_census_gb_par, GVA_per_sector_per_zone_gb, left_on=['WSZONEID'], right_on=['WSZONEID'],how='left')
# wsz_census_gb_par.to_csv('/soge-home/projects/mistral/jamaica-ccri/drought/calibrated_method/input/wsz_census_gb_par_update.csv')

# Water_Supply_Zones = gpd.read_file('potable/Water_Supply_Zone.gpkg')
# Water_Supply_Zones = Water_Supply_Zones.merge(wsz_census_gb_par[['WSZONEID','system_connected_pop_2010','system_connected_pop_2020','system_connected_pop_2030','system_connected_pop_2040','system_connected_pop_2050','GVA (JD/day)','osm_id']], on='WSZONEID') #
# Water_Supply_Zones = gpd.GeoDataFrame(Water_Supply_Zones, crs="EPSG:3448",geometry=Water_Supply_Zones['geometry'])
# Water_Supply_Zones.to_file('Water_Supply_Zone_update.gpkg', driver='GPKG')

# parish_stats = pd.merge(parish_stats,wsz_census_gb[['PARISH','par']],on='PARISH', how='right').drop_duplicates()

# wsz_potable_assets = pd.read_csv('potable/processed/wsz_potable_supply_assets_3.csv')
# wsz_potable_assets['OBJECTID'] = pd.to_numeric(wsz_potable_assets['OBJECTID'])
# wsz_potable_assets['node_id'] = wsz_potable_assets['Type'].astype(str)+'_'+wsz_potable_assets['OBJECTID'].astype(float).astype(str)

# wsz_potable_assets_merge = wsz_potable_assets.merge(wsz_census_gb_par, on='WSZONEID',how='left').drop_duplicates()
# wsz_potable_assets_merge['Type'] = np.where(wsz_potable_assets_merge['LOCATION'].isin(['Mona Dam', 'Hermitage Dam']), 'dam', wsz_potable_assets_merge['Type'])

# source = ['Spring','Treatment Plant','Intake','Production Well','River Source','Filter Plant','Entombment', 'dam']
# junction = ['Relift Station','Pump Station','Booster Station']
# end = ['Storage Tank','Reservoir']
# wsz_potable_assets_merge['node_type'] = np.where(wsz_potable_assets_merge['Type'].isin(source),'source',
                                        # np.where(wsz_potable_assets_merge['Type'].isin(junction),'junction',
                                        # np.where(wsz_potable_assets_merge['Type'].isin(end),'end','nan')))
# wsz_potable_assets_merge_gb = wsz_potable_assets_merge.groupby(['WSZONEID','node_type'])['LOCATION'].count().to_frame().reset_index()
# wsz_potable_assets_merge_gb_end = wsz_potable_assets_merge_gb[wsz_potable_assets_merge_gb['node_type']=='end']
# wsz_potable_assets_merge_gb_end = wsz_potable_assets_merge_gb_end.rename(columns = {'LOCATION':'no. ends'})

# potable_assets = wsz_potable_assets_merge.merge(wsz_potable_assets_merge_gb_end[['WSZONEID','no. ends']], on='WSZONEID',how='left').drop_duplicates()
# potable_assets['frac_pop'] = np.where(potable_assets['Type'].isin(source),1,
                                        # np.where(potable_assets['Type'].isin(junction),1,
                                        # np.where(potable_assets['Type'].isin(end),1/potable_assets['no. ends'],0)))
# potable_assets['asset_connected_pop_2010'] = potable_assets['frac_pop']*potable_assets['system_connected_pop_2010']
# potable_assets['asset_connected_GVA_other_sectors_JD/day'] = potable_assets['frac_pop']*potable_assets['GVA (JD/day)']
# potable_assets['asset_connected_pop_2010'] = potable_assets['asset_connected_pop_2010']/1000
# potable_assets['asset_connected_pop_2010'] = potable_assets['asset_connected_pop_2010'].apply(np.ceil) * 1000

# potable_assets['par'] = potable_assets['WSZONEID'].str[:3]

# # potable_assets = pd.merge(potable_assets, population_growth_zone, left_on='WSZONEID', right_on='WSZONEID')
# potable_assets['asset_connected_pop_2020'] = potable_assets['asset_connected_pop_2010']+potable_assets['asset_connected_pop_2010']*potable_assets['2011_2020_growth']
# potable_assets['asset_connected_pop_2030'] = potable_assets['asset_connected_pop_2010']+potable_assets['asset_connected_pop_2010']*potable_assets['2011_2030_growth']
# potable_assets['asset_connected_pop_2040'] = potable_assets['asset_connected_pop_2010']+potable_assets['asset_connected_pop_2010']*potable_assets['2011_2040_growth']
# potable_assets['asset_connected_pop_2050'] = potable_assets['asset_connected_pop_2010']+potable_assets['asset_connected_pop_2010']*potable_assets['2011_2050_growth']

# potable_assets = pd.merge(potable_assets, parish_stats, left_on='par', right_on='par')

# # # # ### add water scetor GVA and for irrigation 1/5*water and san
# # # # # GVA_per_sector_per_parish = pd.read_csv('potable/processed/GVA_per_zone_per_sector_per_parish.csv') # GVA in 1000
# # # # # GVA_per_sector_per_parish['par'] = GVA_per_sector_per_parish['WSZONEID'].str[:3]
# # # # # GVA_per_sector_per_parish['building_type_nat_GVA'] = np.where(GVA_per_sector_per_parish['building_type'].isin(['Resort','Recreation']),'hotels/condominium',
                                                  # # # # # np.where(GVA_per_sector_per_parish['building_type'].isin(['Commercial']),'commerical',
                                                  # # # # # np.where(GVA_per_sector_per_parish['building_type'].isin(['Institutional']),'gov/school','nan')))

# # # # # GVA_per_sector_per_parish_gb_par = GVA_per_sector_per_parish.groupby(['par','building_type_nat_GVA'])['osm_id'].count().to_frame().reset_index()
# # # # # GVA_per_sector_per_parish_gb_par = GVA_per_sector_per_parish_gb_par.rename(columns = {'osm_id':'par_count'})

# # # # # GVA_per_sector_per_zone_gb = GVA_per_sector_per_parish.groupby(['WSZONEID','building_type_nat_GVA'])['osm_id'].apply(list).to_frame().reset_index()
# # # # # GVA_per_sector_per_zone_gb = GVA_per_sector_per_zone_gb.rename(columns = {'osm_id':'list_osm_id'})
# # # # # GVA_per_sector_per_zone_gb['count'] = GVA_per_sector_per_parish.groupby(['par','WSZONEID','building_type_nat_GVA'])['osm_id'].count().to_frame().reset_index()['osm_id']
# # # # # GVA_per_sector_per_zone_gb = GVA_per_sector_per_zone_gb.merge(GVA_per_sector_per_parish_gb_par, on=['par','building_type_nat_GVA'])
# # # # # GVA_per_sector_per_zone_gb['frac_GVA_sector_zone_per_par'] = GVA_per_sector_per_zone_gb['count']/GVA_per_sector_per_zone_gb['par_count']
# # # # # GVA_per_sector_per_parish_gb_condominium = GVA_per_sector_per_zone_gb[GVA_per_sector_per_zone_gb['building_type_nat_GVA']=='hotels/condominium']
# # # # # GVA_per_sector_per_parish_gb_condominium = GVA_per_sector_per_parish_gb_condominium.rename(columns = {'list_osm_id':'list_condominium','frac_GVA_sector_zone_per_par':'frac_condominium'})
# # # # # GVA_per_sector_per_parish_gb_commercial = GVA_per_sector_per_zone_gb[GVA_per_sector_per_zone_gb['building_type_nat_GVA']=='commerical']
# # # # # GVA_per_sector_per_parish_gb_commercial = GVA_per_sector_per_parish_gb_commercial.rename(columns = {'list_osm_id':'list_commercial','frac_GVA_sector_zone_per_par':'frac_commercial'})
# # # # # GVA_per_sector_per_parish_gb_government = GVA_per_sector_per_zone_gb[GVA_per_sector_per_zone_gb['building_type_nat_GVA']=='gov/school']
# # # # # GVA_per_sector_per_parish_gb_government = GVA_per_sector_per_parish_gb_government.rename(columns = {'list_osm_id':'list_government','frac_GVA_sector_zone_per_par':'frac_government'})

# # # # # potable_assets = potable_assets.merge(GVA_per_sector_per_parish_gb_condominium, on=['par','WSZONEID'])
# # # # # potable_assets = potable_assets.merge(GVA_per_sector_per_parish_gb_commercial, on=['par','WSZONEID'])
# # # # # potable_assets = potable_assets.merge(GVA_per_sector_per_parish_gb_government, on=['par','WSZONEID'])

# # # # # potable_assets['asset_GVA - hotels/condominium'] = potable_assets['frac_pop'] * potable_assets['GVA - hotels/condominium'] * potable_assets['frac_condominium'] # frac GVA - hotels/condominium per system within parish
# # # # # potable_assets['asset_GVA - commerical'] = potable_assets['frac_pop'] * potable_assets['GVA - commerical'] * potable_assets['frac_commercial']
# # # # # potable_assets['asset_GVA - gov/school'] = potable_assets['frac_pop'] * potable_assets['GVA - gov/school'] * potable_assets['frac_government']
# # # # # potable_assets['asset_GVA - total'] = potable_assets['frac_pop'] * potable_assets['GVA']
# # # # # potable_assets.to_csv('potable_assets.csv')

# # # # # check = potable_assets.groupby('WSZONEID')['asset_connected_pop_2020'].max().to_frame().reset_index()
# # # # # print(check['asset_connected_pop_2020'].sum())
# # # # # print(oliv)
# # # # # potable_assets_gb['asset_GVA - hotels/condominium'] = potable_assets.groupby('WSZONEID')['asset_GVA - hotels/condominium'].max().to_frame().reset_index()['asset_GVA - hotels/condominium']

# potable_facilities_NWC = gpd.read_file('/soge-home/projects/mistral/jamaica-ccri/processed_data/networks/water/potable_facilities_NWC.gpkg')
# potable_facilities_NWC = potable_facilities_NWC.drop(['asset_pop_new_2020'], axis=1)
# potable_facilities_NWC_pop = potable_facilities_NWC.merge(potable_assets[['node_id','node_type','X','Y','WSZONEID','asset_connected_pop_2010','asset_connected_pop_2020','asset_connected_pop_2030','asset_connected_pop_2040','asset_connected_pop_2050','asset_connected_GVA_other_sectors_JD/day','osm_id']], on='node_id', how='left').drop_duplicates()

# potable_facilities_NWC_pop['W_GDP'] = 17648*(1e6/365)*0.8 * (potable_facilities_NWC_pop['asset_connected_pop_2020']/potable_facilities_NWC_pop['asset_connected_pop_2020'].sum())

# potable_facilities_NWC_economic = potable_facilities_NWC_pop[(potable_facilities_NWC_pop['node_type'].isin(['source','junction']))&(potable_facilities_NWC_pop['asset_connected_pop_2020']>0)]
# potable_facilities_NWC_economic_m = potable_facilities_NWC_economic[['node_id','WSZONEID']].merge(GVA_per_sector_per_parish[['osm_id','A_GDP', 'B_GDP', 'C_GDP', 'D_GDP', 'E_GDP', 'F_GDP', 'G_GDP', 'H_GDP',
       # 'I_GDP', 'J_GDP', 'K_GDP', 'L_GDP', 'M_GDP', 'N_GDP', 'O_GDP',
       # 'total_GDP', 'GDP_unit','WSZONEID']], on='WSZONEID', how='left')
# potable_facilities_NWC_economic_m = potable_facilities_NWC_economic_m[['node_id','osm_id','A_GDP', 'B_GDP', 'C_GDP', 'D_GDP', 'E_GDP', 'F_GDP', 'G_GDP', 'H_GDP',
       # 'I_GDP', 'J_GDP', 'K_GDP', 'L_GDP', 'M_GDP', 'N_GDP', 'O_GDP',
       # 'total_GDP', 'GDP_unit','WSZONEID']]
# potable_facilities_NWC_economic_m.to_csv('potable_facilities_buildings_economic_activity_mapping.csv') # W_GDP

# potable_facilities_NWC_economic_m_2 = potable_facilities_NWC_pop[['node_id','W_GDP','WSZONEID']].merge(GVA_per_sector_per_parish[['osm_id','A_GDP', 'B_GDP', 'C_GDP', 'D_GDP', 'E_GDP', 'F_GDP', 'G_GDP', 'H_GDP',
       # 'I_GDP', 'J_GDP', 'K_GDP', 'L_GDP', 'M_GDP', 'N_GDP', 'O_GDP',
       # 'total_GDP', 'GDP_unit','WSZONEID']], on='WSZONEID', how='left')
# potable_facilities_NWC_economic_m_2_gb = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['W_GDP'].mean().to_frame().reset_index()
# potable_facilities_NWC_economic_m_2_gb['A_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['A_GDP'].sum().to_frame().reset_index()['A_GDP']
# potable_facilities_NWC_economic_m_2_gb['B_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['B_GDP'].sum().to_frame().reset_index()['B_GDP']
# potable_facilities_NWC_economic_m_2_gb['C_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['C_GDP'].sum().to_frame().reset_index()['C_GDP']
# potable_facilities_NWC_economic_m_2_gb['D_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['D_GDP'].sum().to_frame().reset_index()['D_GDP']
# potable_facilities_NWC_economic_m_2_gb['E_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['E_GDP'].sum().to_frame().reset_index()['E_GDP']
# potable_facilities_NWC_economic_m_2_gb['F_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['F_GDP'].sum().to_frame().reset_index()['F_GDP']
# potable_facilities_NWC_economic_m_2_gb['G_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['G_GDP'].sum().to_frame().reset_index()['G_GDP']
# potable_facilities_NWC_economic_m_2_gb['H_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['H_GDP'].sum().to_frame().reset_index()['H_GDP']
# potable_facilities_NWC_economic_m_2_gb['I_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['I_GDP'].sum().to_frame().reset_index()['I_GDP']
# potable_facilities_NWC_economic_m_2_gb['J_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['J_GDP'].sum().to_frame().reset_index()['J_GDP']
# potable_facilities_NWC_economic_m_2_gb['K_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['K_GDP'].sum().to_frame().reset_index()['K_GDP']
# potable_facilities_NWC_economic_m_2_gb['L_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['L_GDP'].sum().to_frame().reset_index()['L_GDP']
# potable_facilities_NWC_economic_m_2_gb['M_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['M_GDP'].sum().to_frame().reset_index()['M_GDP']
# potable_facilities_NWC_economic_m_2_gb['N_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['N_GDP'].sum().to_frame().reset_index()['N_GDP']
# potable_facilities_NWC_economic_m_2_gb['O_GDP'] = potable_facilities_NWC_economic_m_2.groupby(['node_id','WSZONEID'])['O_GDP'].sum().to_frame().reset_index()['O_GDP']
    # #'A_GDP', 'B_GDP', 'C_GDP', 'D_GDP', 'E_GDP', 'F_GDP', 'G_GDP', 'H_GDP',
    # #'I_GDP', 'J_GDP', 'K_GDP', 'L_GDP', 'M_GDP', 'N_GDP', 'O_GDP'   
# print(potable_facilities_NWC_economic_m_2_gb)
# potable_facilities_NWC_economic_m_2_gb.to_csv('potable_facilities_dependent_economic_activity.csv')
# ## print(oliv)
# potable_facilities_NWC_pop['cost ($J) - lower bound'] = np.where(potable_facilities_NWC_pop['asset_type']=='well',5000000,potable_facilities_NWC_pop['cost ($J) - lower bound']) # same for irrigation
# potable_facilities_NWC_pop['cost ($J) - upper bound'] = np.where(potable_facilities_NWC_pop['asset_type']=='well',20000000,potable_facilities_NWC_pop['cost ($J) - upper bound']) # same for irrigation
# potable_facilities_NWC_pop['min_damage_cost'] = np.where(potable_facilities_NWC_pop['asset_type']=='well',5000000,potable_facilities_NWC_pop['min_damage_cost']) # same for irrigation
# potable_facilities_NWC_pop['max_damage_cost'] = np.where(potable_facilities_NWC_pop['asset_type']=='well',20000000,potable_facilities_NWC_pop['max_damage_cost']) # same for irrigation
# # potable_facilities_NWC = potable_facilities_NWC.drop(['asset_pop_new_2010','asset_pop_new_2020','asset_pop_new_2030'], axis=1)
# # potable_facilities_NWC_pop = potable_facilities_NWC.merge(potable_assets[['node_id','WSZONEID','asset_connected_pop_2010']], on='node_id', how='left').drop_duplicates()
# # potable_facilities_NWC_pop['asset_connected_pop_2022'] = potable_facilities_NWC_pop['asset_connected_pop_2022'].apply(roundup)

# potable_facilities_NWC_pop.to_file('potable_facilities_NWC_pop_update.gpkg', driver='GPKG')

# # # # # potable_assets = potable_facilities_NWC_pop
# # # # # pipelines_centroid = pd.read_csv('potable/processed/pipelines_centroid.csv')
# # # # # pipelines_centroid['OBJECTID'] = pd.to_numeric(pipelines_centroid['OBJECTID'])
# # # # # pipelines_centroid['edge_id'] = 'potable_'+pipelines_centroid['OBJECTID'].astype(float).astype(str)
# # # # # for it,row in pipelines_centroid.iterrows(): #canal/pipe centroid
    # # # # # print(row)
    # # # # # pipeline_id = row.OBJECTID #OBJECTID
    # # # # # pipelines_Y = row.Y
    # # # # # pipelines_X = row.X
    # # # # # potable_assets['lat_diff'] = abs(potable_assets['Y'] - pipelines_Y)
    # # # # # potable_assets['lon_diff'] = abs(potable_assets['X'] - pipelines_X)
    # # # # # potable_assets['sqrt'] = np.sqrt(np.square(potable_assets['lat_diff'])+np.square(potable_assets['lat_diff']))
    # # # # # closest_wrz = potable_assets[potable_assets['sqrt']==potable_assets['sqrt'].min()]['WSZONEID'].values[0]
    # # # # # closest_node = potable_assets[potable_assets['sqrt']==potable_assets['sqrt'].min()]['node_id'].values[0]
    # # # # # asset_connected_pop_2010 = potable_assets[potable_assets['sqrt']==potable_assets['sqrt'].min()]['asset_connected_pop_2010'].values[0]
    # # # # # asset_connected_pop_2020 = potable_assets[potable_assets['sqrt']==potable_assets['sqrt'].min()]['asset_connected_pop_2020'].values[0]
    # # # # # asset_connected_pop_2030 = potable_assets[potable_assets['sqrt']==potable_assets['sqrt'].min()]['asset_connected_pop_2030'].values[0]
    # # # # # asset_connected_pop_2040 = potable_assets[potable_assets['sqrt']==potable_assets['sqrt'].min()]['asset_connected_pop_2040'].values[0]
    # # # # # asset_connected_pop_2050 = potable_assets[potable_assets['sqrt']==potable_assets['sqrt'].min()]['asset_connected_pop_2050'].values[0]
    # # # # # #asset_connected_GVA_other_sectors = potable_assets[potable_assets['sqrt']==potable_assets['sqrt'].min()]['asset_connected_GVA_other_sectors_JD/day'].values[0]
    # # # # # #asset_connected_GVA_water_sector = potable_assets[potable_assets['sqrt']==potable_assets['sqrt'].min()]['asset_connected_GVA_water_sector_JD/day'].values[0]
    # # # # # pipelines_centroid.loc[it,'WSZONEID'] = closest_wrz
    # # # # # pipelines_centroid.loc[it,'node_id'] = closest_node
    # # # # # pipelines_centroid.loc[it,'asset_connected_pop_2010'] = asset_connected_pop_2010
    # # # # # pipelines_centroid.loc[it,'asset_connected_pop_2020'] = asset_connected_pop_2020
    # # # # # pipelines_centroid.loc[it,'asset_connected_pop_2030'] = asset_connected_pop_2030
    # # # # # pipelines_centroid.loc[it,'asset_connected_pop_2040'] = asset_connected_pop_2040
    # # # # # pipelines_centroid.loc[it,'asset_connected_pop_2050'] = asset_connected_pop_2050
    # # # # # #pipelines_centroid.loc[it,'asset_connected_GVA_other_sectors'] = asset_connected_GVA_other_sectors
    # # # # # #pipelines_centroid.loc[it,'asset_connected_GVA_water_sector'] = asset_connected_GVA_water_sector
    
# pipelines_centroid = gpd.read_file('potable/processed/pipelines_centroid_filtered.gpkg')

# pipelines_NWC = gpd.read_file('/soge-home/projects/mistral/jamaica-ccri/processed_data/networks/water/pipelines_NWC.gpkg')
# pipelines_NWC_pop = pipelines_NWC.merge(pipelines_centroid[['edge_id','node_id']], on='edge_id', how='left').drop_duplicates()
# # # pipelines_NWC_pop = pipelines_NWC_pop.dropna(subset=['asset_connected_pop_2010'])
# # # pipelines_NWC_pop['asset_connected_pop_2010'] = pipelines_NWC_pop['asset_connected_pop_2010'].apply(roundup)
# pipelines_NWC_pop.to_file('pipelines_NWC_edges_pop_update.gpkg', driver='GPKG')
# pipelines_NWC_pop = gpd.read_file('pipelines_NWC_edges_pop_update.gpkg')
# potable_pipelines_NWC_economic_m = pipelines_NWC_pop[['edge_id','node_id']].merge(potable_facilities_NWC_economic_m[['node_id','osm_id','A_GDP', 'B_GDP', 'C_GDP', 'D_GDP', 'E_GDP', 'F_GDP', 'G_GDP', 'H_GDP',
       # 'I_GDP', 'J_GDP', 'K_GDP', 'L_GDP', 'M_GDP', 'N_GDP', 'O_GDP',
       # 'total_GDP', 'GDP_unit','WSZONEID']], on='node_id', how='left')
# potable_pipelines_NWC_economic_m[['edge_id','osm_id','A_GDP', 'B_GDP', 'C_GDP', 'D_GDP', 'E_GDP', 'F_GDP', 'G_GDP', 'H_GDP',
       # 'I_GDP', 'J_GDP', 'K_GDP', 'L_GDP', 'M_GDP', 'N_GDP', 'O_GDP',
       # 'total_GDP', 'GDP_unit']].to_csv('potable_pipelines_buildings_economic_activity_mapping.csv') # W_GDP
# potable_pipelines_NWC_economic_m = pipelines_NWC_pop[['edge_id','node_id']].merge(potable_facilities_NWC_economic_m_2_gb[['node_id','W_GDP','A_GDP', 'B_GDP', 'C_GDP', 'D_GDP', 'E_GDP', 'F_GDP', 'G_GDP', 'H_GDP',
       # 'I_GDP', 'J_GDP', 'K_GDP', 'L_GDP', 'M_GDP', 'N_GDP', 'O_GDP']], on='node_id', how='left').to_csv('potable_pipelines_dependent_economic_activity.csv')
# print(potable_pipelines_NWC_economic_m)


########################################################################
########################################################################
# disruption analysis potable
########################################################################
########################################################################

# damaged_nodes_df = pd.DataFrame([])
# merged_nodes = pd.merge(potable_assets,damaged_nodes_df,on='node_id')
# pop_disrupted_nodes = merged_nodes.groupby('WSZONEID')['asset_connected_pop_2010'].max().reset_index().to_frame() ### sink assets are storage tanks that are assumed not to be vulnerable
# ## same for GVA

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

# agriculture = pd.read_csv('irrigation/processed/agriculture_gdp.csv') #merge with schemes
# print(agriculture,agriculture.columns)
# agriculture['GDP_persqm'] = pd.to_numeric(agriculture['GDP_persqm'], errors='coerce')
# agriculture['area_m2'] = pd.to_numeric(agriculture['area_m2'], errors='coerce')
# agriculture['GVA (JD/day)'] = agriculture['GDP_persqm']*agriculture['area_m2']
# agriculture['DXF_TEXT'] = np.where(agriculture['District']=='MID - CLARENDON', 'MID CLARENDON',
                          # np.where(agriculture['District']=='RIO - COBRE', 'Rio Cobre', agriculture['District']))
# agriculture_gb = agriculture.groupby('DXF_TEXT')['GVA (JD/day)'].sum().to_frame().reset_index()
# agriculture_gb['GDP_persqm'] = agriculture.groupby('DXF_TEXT')['GDP_persqm'].mean().to_frame().reset_index()['GDP_persqm']
# agriculture_gb['economic_interection_area'] = agriculture.groupby('DXF_TEXT')['area_m2'].sum().to_frame().reset_index()['area_m2']

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
agriculture['intersect_area_sqm'] = agriculture['area'] #* 1e-3 * 1e-3
agriculture['GDP_persqm'] = pd.to_numeric(agriculture['GDP_persqm'], errors='coerce')
agriculture['GVA (JD/day)'] = agriculture['GDP_persqm']*agriculture['intersect_area_sqm']
agriculture['DXF_TEXT'] = np.where(agriculture['District']=='MID - CLARENDON', 'MID CLARENDON',
                          np.where(agriculture['District']=='RIO - COBRE', 'Rio Cobre', agriculture['District']))
agriculture_gb = agriculture.groupby('DXF_TEXT')['GVA (JD/day)'].sum().to_frame().reset_index()
agriculture_gb['GDP_persqm'] = agriculture.groupby('DXF_TEXT')['GDP_persqm'].mean().to_frame().reset_index()['GDP_persqm']
agriculture_gb['economic_interection_area'] = agriculture.groupby('DXF_TEXT')['intersect_area_sqm'].sum().to_frame().reset_index()['intersect_area_sqm']

schemes = NIC_IRRIGATION_SCHEMES.merge(agriculture_gb, on='DXF_TEXT')
# schemes['GVA (JD/day)'] = schemes['GVA (JD/day)']*(schemes['scheme_area_sqm']/schemes['economic_interection_area'])
schemes_gb = schemes.groupby(['DXF_TEXT'])['GVA (JD/day)'].sum().to_frame().reset_index()
print(schemes_gb, NIC_IRRIGATION_SCHEMES)
NIC_IRRIGATION_SCHEMES = NIC_IRRIGATION_SCHEMES.merge(schemes_gb, on='DXF_TEXT')
NIC_IRRIGATION_SCHEMES = gpd.GeoDataFrame(NIC_IRRIGATION_SCHEMES,crs="EPSG:3448",geometry=NIC_IRRIGATION_SCHEMES['geometry']) #.merge(schemes_gb, on='DXF_TEXT')
print(NIC_IRRIGATION_SCHEMES)
NIC_IRRIGATION_SCHEMES.to_file('NIC_IRRIGATION_SCHEMES_update.gpkg', driver='GPKG')

wells = pd.read_csv('irrigation/wells.csv')
wells['OBJECTID'] = np.linspace(0,wells.shape[0],wells.shape[0])
wells['node_id'] = 'well_'+wells['OBJECTID'].astype(float).astype(str)
wells_node_id = wells.to_csv('irrigation/wells_node_id.csv')

# scheme_centroids = pd.read_csv('irrigation/processed/schemes_centroids.csv')
# scheme_centroids = scheme_centroids.merge(schemes_gb, on='DXF_TEXT')
# pipelines_centroids = pd.read_csv('irrigation/processed/pipe_centroids.csv')
# pipelines_centroids['edge_id'] = 'pipeline_'+pipelines_centroids['OBJECTID'].astype(float).astype(str) 

canals_centroids = pd.read_csv('irrigation/processed/canal_centroids.csv')

wells['DXF_TEXT'] = np.where(wells['DISTRICT']=='MID - CLARENDON', 'MID CLARENDON',
                    np.where(wells['DISTRICT']=='RIO - COBRE', 'Rio Cobre', wells['DISTRICT']))

wells_NIC_economic_m = wells[['node_id','DXF_TEXT']].merge(agriculture[['fid','land_id','forest_id','GDP_persqm', 'intersect_area_sqm', 'GDP_unit', 'DXF_TEXT']], on='DXF_TEXT', how='left')
wells_NIC_economic_m['GVA (JD/day)'] = wells_NIC_economic_m['GDP_persqm']*wells_NIC_economic_m['intersect_area_sqm']

wells_NIC_economic_m[['node_id','fid','land_id','forest_id','GDP_persqm', 'intersect_area_sqm', 'GDP_unit', 'GVA (JD/day)', 'DXF_TEXT']].to_csv('irrigation_nodes_agriculture_economic_activity_mapping.csv') # W_GDP
print(wells_NIC_economic_m)
wells_NIC_economic_m_2 = wells[['node_id','DXF_TEXT']].merge(agriculture[['fid','land_id','forest_id','GDP_persqm', 'intersect_area_sqm', 'GDP_unit', 'DXF_TEXT']], on='DXF_TEXT', how='left')
wells_NIC_economic_m_2['GVA (JD/day)'] = wells_NIC_economic_m_2['GDP_persqm']*wells_NIC_economic_m_2['intersect_area_sqm']

wells_NIC_economic_m_3 = wells_NIC_economic_m_2.groupby('node_id')['GVA (JD/day)'].sum().to_frame().reset_index()  
print(wells_NIC_economic_m_3)
wells_NIC_economic_m_3.to_csv('irrigation_nodes_dependent_economic_activity.csv')
# print(oliv)

wells_merge = pd.merge(wells,NIC_IRRIGATION_SCHEMES,on='DXF_TEXT',how='left').drop_duplicates()
wells_merge_gb = wells_merge.groupby(['DXF_TEXT','scheme_area_sqm','GVA (JD/day)'])['NAME'].count().to_frame().reset_index()
wells_merge_gb = wells_merge_gb.rename(columns={'NAME':'no. wells'})
wells_merge_2 = wells_merge.merge(wells_merge_gb, on=['DXF_TEXT','scheme_area_sqm','GVA (JD/day)'],how='left').drop_duplicates()
wells_merge_2['sqm_per_node'] = wells_merge_2['scheme_area_sqm']/wells_merge_2['no. wells']
wells_merge_2['GDP/day_per_node'] = wells_merge_2['GVA (JD/day)']/wells_merge_2['no. wells']
#wells_merge_2['acres_per_node'] = wells_merge_2['acres_per_node'].apply(roundup)
print(wells_merge_2.columns)
irrigation_assets_NIC = gpd.read_file('/soge-home/projects/mistral/jamaica-ccri/processed_data/networks/water/irrigation_assets_NIC.gpkg', layer ='nodes')
irrigation_assets_NIC['node_id'] = 'well_'+irrigation_assets_NIC['OBJECTID'].astype(float).astype(str)
irrigation_assets_NIC['cost ($J) - lower bound'] = np.where(irrigation_assets_NIC['asset_type']=='well',5000000,irrigation_assets_NIC['cost ($J) - lower bound']) # same for irrigation
irrigation_assets_NIC['cost ($J) - upper bound'] = np.where(irrigation_assets_NIC['asset_type']=='well',20000000,irrigation_assets_NIC['cost ($J) - upper bound']) # same for irrigation
irrigation_assets_NIC['min_damage_cost'] = np.where(irrigation_assets_NIC['asset_type']=='well',5000000,irrigation_assets_NIC['min_damage_cost']) # same for irrigation
irrigation_assets_NIC['max_damage_cost'] = np.where(irrigation_assets_NIC['asset_type']=='well',20000000,irrigation_assets_NIC['max_damage_cost']) # same for irrigation
irrigation_assets_NIC_are = irrigation_assets_NIC.merge(wells_merge_2[['node_id','DXF_TEXT', 'sqm_per_node','GDP/day_per_node']], on='node_id', how='left').drop_duplicates()
irrigation_assets_NIC_are.to_file('irrigation_assets_NIC_are_nodes_update.gpkg', driver='GPKG')

# # # # # # for it,row in pipelines_centroids.iterrows():
    # # # # # # pipeline_id = row.OBJECTID 
    # # # # # # pipelines_Y = row.Y
    # # # # # # pipelines_X = row.X
    # # # # # # pipelines_area = row.Size
    # # # # # # scheme_centroids['lat_diff'] = abs(scheme_centroids['Y'] - pipelines_Y)
    # # # # # # scheme_centroids['lon_diff'] = abs(scheme_centroids['X'] - pipelines_X)
    # # # # # # scheme_centroids['sqrt'] = np.sqrt(np.square(scheme_centroids['lat_diff'])+np.square(scheme_centroids['lat_diff']))
    # # # # # # closest_scheme = scheme_centroids[scheme_centroids['sqrt']==scheme_centroids['sqrt'].min()]['DXF_TEXT'].values[0]
    # # # # # # closest_node = scheme_centroids[scheme_centroids['sqrt']==scheme_centroids['sqrt'].min()]['node_id'].values[0]
    # # # # # # pipelines_centroids.loc[it,'DXF_TEXT'] = closest_scheme
    # # # # # # pipelines_centroids.loc[it,'node_id'] = closest_node
    # # # # # # scheme_area_sqm = scheme_centroids[scheme_centroids['sqrt']==scheme_centroids['sqrt'].min()]['scheme_area_sqm'].values[0]
    # # # # # # pipelines_centroids.loc[it,'scheme_area_sqm'] = scheme_area_sqm
    # # # # # # GDP_persqm = scheme_centroids[scheme_centroids['sqrt']==scheme_centroids['sqrt'].min()]['GDP_persqm'].values[0]
    # # # # # # pipelines_centroids.loc[it,'GDP_persqm'] = GDP_persqm

pipelines_centroids = gpd.read_file('irrigation/processed/pipe_centroids_output_from_qgis_2.gpkg') ## read from qgis
pipelines_centroids['edge_id'] = 'pipeline_'+pipelines_centroids['OBJECTID'].astype(float).astype(str) 
pipelines_centroids = pipelines_centroids.merge(irrigation_assets_NIC_are[['DXF_TEXT','node_id']], on='node_id')
print(pipelines_centroids.columns)
# print(oliv)
# ## irrigation_assets_NIC_are merge to get DX T tet

# pipelines_centroids_merge = pd.merge(pipelines_centroids,NIC_IRRIGATION_SCHEMES[['DXF_TEXT','scheme_area_sqm','GVA (JD/day)']],on=['DXF_TEXT'],how='left').drop_duplicates()
# pipelines_centroids_merge_gb = pipelines_centroids_merge.groupby(['DXF_TEXT','scheme_area_sqm'])['Size'].sum().to_frame().reset_index()
# pipelines_centroids_merge_gb = pipelines_centroids_merge_gb.rename(columns={'Size':'sum of pipe size'})
# pipelines_centroids_merge_2 = pd.merge(pipelines_centroids_merge, pipelines_centroids_merge_gb, on=['DXF_TEXT','scheme_area_sqm'], how='left').drop_duplicates()
# pipelines_centroids_merge_2['sqm_per_edge'] = pipelines_centroids_merge_2['scheme_area_sqm']*pipelines_centroids_merge_2['Size']/pipelines_centroids_merge_2['sum of pipe size']
# pipelines_centroids_merge_2['GDP/day_per_edge'] = pipelines_centroids_merge_2['GVA (JD/day)']*pipelines_centroids_merge_2['Size']/pipelines_centroids_merge_2['sum of pipe size']
# # pipelines_centroids_merge_2['GDP/day_per_edge'] = pipelines_centroids_merge_2['sqm_per_edge']*pipelines_centroids_merge_2['GDP_persqm']

# for it,row in canals_centroids.iterrows(): #canal/pipe centroid
    # pipeline_id = row.OBJECTID #
    # pipelines_Y = row.Y
    # pipelines_X = row.X
    # scheme_centroids['lat_diff'] = abs(scheme_centroids['Y'] - pipelines_Y)
    # scheme_centroids['lon_diff'] = abs(scheme_centroids['X'] - pipelines_X)
    # scheme_centroids['sqrt'] = np.sqrt(np.square(scheme_centroids['lat_diff'])+np.square(scheme_centroids['lat_diff']))
    # closest_scheme = scheme_centroids[scheme_centroids['sqrt']==scheme_centroids['sqrt'].min()]['DXF_TEXT'].values[0]
    # closest_node = scheme_centroids[scheme_centroids['sqrt']==scheme_centroids['sqrt'].min()]['node_id'].values[0]
    # canals_centroids.loc[it,'DXF_TEXT'] = closest_scheme
    # pipelines_centroids.loc[it,'node_id'] = closest_node
    # scheme_area_sqm = scheme_centroids[scheme_centroids['sqrt']==scheme_centroids['sqrt'].min()]['scheme_area_sqm'].values[0]
    # canals_centroids.loc[it,'scheme_area_sqm'] = scheme_area_sqm
    # GDP_persqm = scheme_centroids[scheme_centroids['sqrt']==scheme_centroids['sqrt'].min()]['GDP_persqm'].values[0]
    # canals_centroids.loc[it,'GDP_persqm'] = GDP_persqm
    
canals_centroids = gpd.read_file('irrigation/processed/canal_centroids_output_from_qgis_2.gpkg') ## read from qgis
canals_centroids['edge_id'] = 'pipeline_'+canals_centroids['OBJECTID'].astype(float).astype(str) 
canals_centroids = canals_centroids.merge(irrigation_assets_NIC_are[['DXF_TEXT','node_id']], on='node_id')

# canals_centroids['DXF_TEXT'] = np.where(canals_centroids['Parish']=='Clarendon','MID CLARENDON',canals_centroids['DXF_TEXT'])
# canals_centroids_merge = pd.merge(canals_centroids,NIC_IRRIGATION_SCHEMES,on=['DXF_TEXT'],how='left').drop_duplicates()
# canals_centroids_merge_gb = canals_centroids_merge.groupby(['DXF_TEXT','scheme_area_sqm'])['LENGTH'].sum().to_frame().reset_index()
# canals_centroids_merge_gb = canals_centroids_merge_gb.rename(columns={'LENGTH':'sum of canal length'})
# canals_centroids_merge_2 = canals_centroids_merge.merge(canals_centroids_merge_gb, on=['DXF_TEXT','scheme_area_sqm'],how='left').drop_duplicates()
# canals_centroids_merge_2['sqm_per_edge'] = canals_centroids_merge_2['scheme_area_sqm']*canals_centroids_merge_2['LENGTH']/canals_centroids_merge_2['sum of canal length']
# canals_centroids_merge_2['GDP/day_per_edge'] = canals_centroids_merge_2['GVA (JD/day)']*canals_centroids_merge_2['LENGTH']/canals_centroids_merge_2['sum of canal length']
# print(canals_centroids_merge_2.columns)
# canals_centroids_merge_2['GDP/day_per_edge'] = canals_centroids_merge_2['sqm_per_edge']*canals_centroids_merge_2['GDP_persqm_y']

# canals_centroids_merge_2['sqm_per_edge'] = np.where(canals_centroids_merge_2['DXF_TEXT']=='MID CLARENDON',canals_centroids_merge_2['AREA_SERVE'],canals_centroids_merge_2['sqm_per_edge'])
edges_centroids = pd.concat([pipelines_centroids,canals_centroids])
#edges_centroids = edges_centroids.dropna(subset=['acres_per_edge'])
#edges_centroids['acres_per_edge'] = edges_centroids['acres_per_edge'].apply(roundup)

irrigation_assets_NIC_edges = gpd.read_file('irrigation/irrigation_edges.gpkg')
irrigation_assets_NIC_edges['edge_id'] = 'pipeline_'+irrigation_assets_NIC_edges['OBJECTID'].astype(str)
irrigation_assets_NIC_edges['edge_id'] = irrigation_assets_NIC_edges['edge_id'].astype(str)
edges_centroids['edge_id'] = edges_centroids['edge_id'].astype(str)
#['OBJECTID_x', 'asset_type', 'asset_type_cost_data',
       # 'asset_type_flood_damage', 'asset_type_hurricane_damage', 'asset',
       # 'curve', 'cost ($J) - lower bound', 'cost ($J) - upper bound', 'source',
       # 'comment', 'cost_unit', 'min_damage_cost', 'max_damage_cost', 'edge_id',
       # 'geometry_x', 'X', 'Y', 'OBJECTID_y', 'ID', 'LENGTH', 'Size', 'type',
       # 'status_1', 'REGION', 'Shape_Leng', 'field_1', 'NAME', 'SERIAL_NUM',
       # 'METER_TYP', 'DATE_', 'DISTRICT', 'LOCATION', 'X_2', 'Y_2', 'CODE',
       # 'WRA_ID', 'LICENCE_NO', 'ABSTRACT_R', 'NIC_NAME', 'NIC_NAME_1', 'DATE',
       # 'COMMENNTS', 'TYPE_2', 'OBJECTID_2', 'node_id', 'n', 'distance',
       # 'feature_x', 'feature_y', 'nearest_x', 'nearest_y', 'geometry_y',
       # 'DXF_TEXT', 'scheme_area_sqm', 'GVA (JD/day)', 'sum of pipe size',
       # 'sqm_per_edge', 'GDP/day_per_edge', 'FNODE_', 'TNODE_', 'LPOLY_',
       # 'RPOLY_', 'CANAL_', 'DXF_LAYER', 'DXF_COLOR', 'DXF_THICKN', 'DXF_TYPE',
       # 'DXF_ELEVAT', 'DXF_HANDLE', 'DXF_CURVE', 'DIST_NAME', 'CANAL_NAME',
       # 'AREA_ZONE', 'LENGTH__M_', 'LENGTH__FT', 'MATERIAL', 'DESIGN',
       # 'CONDITION', 'TYPE', 'AREA_SERVE', 'CONTR_VOLU', 'MAJOR_CROP',
       # 'WATER_SOUR', 'LENG', 'KM', 'CLASS', 'Parish', 'DISTRICT_2', 'CODE_2',
       # 'sum of canal length'],
irrigation_assets_NIC_are_egdes = irrigation_assets_NIC_edges.merge(edges_centroids[['edge_id','node_id']], on='edge_id', how='left').drop_duplicates()
print(irrigation_assets_NIC_are_egdes.columns)
irrigation_assets_NIC_are_egdes = gpd.GeoDataFrame(irrigation_assets_NIC_are_egdes,crs="EPSG:3448",geometry=irrigation_assets_NIC_are_egdes['geometry'])
irrigation_assets_NIC_are_egdes.to_file('irrigation_assets_NIC_are_edges_update.gpkg', driver='GPKG') #.drop_duplicates()

irrigation_edges_NIC_economic_m = edges_centroids[['edge_id','node_id']].merge(wells_NIC_economic_m[['node_id','fid','land_id','forest_id', 'GVA (JD/day)', 'DXF_TEXT']], on='node_id', how='left')
irrigation_edges_NIC_economic_m.to_csv('irrigation_edges_agriculture_economic_activity_mapping.csv') # W_GDP
irrigation_edges_NIC_economic_m = edges_centroids[['edge_id','node_id']].merge(wells_NIC_economic_m_3, on='node_id', how='left')
irrigation_edges_NIC_economic_m.to_csv('irrigation_edges_dependent_economic_activity.csv')
print(irrigation_edges_NIC_economic_m)
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
    
