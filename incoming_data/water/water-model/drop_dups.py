import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd

pipelines_NWC_edges_pop_update = gpd.read_file('potable_facilities_NWC_pop_update.gpkg')
pipelines_NWC_edges_pop_update = pipelines_NWC_edges_pop_update.drop_duplicates(subset=['node_id','WSZONEID'])
pipelines_NWC_edges_pop_update_gb = pipelines_NWC_edges_pop_update.groupby('node_id')['asset_connected_pop_2020'].sum().to_frame().reset_index()
pipelines_NWC_edges_pop_update_gb['asset_connected_pop_2030'] = pipelines_NWC_edges_pop_update.groupby('node_id')['asset_connected_pop_2030'].sum().to_frame().reset_index()['asset_connected_pop_2030']
pipelines_NWC_edges_pop_update_gb['asset_connected_pop_2040'] = pipelines_NWC_edges_pop_update.groupby('node_id')['asset_connected_pop_2040'].sum().to_frame().reset_index()['asset_connected_pop_2040']
pipelines_NWC_edges_pop_update_gb['asset_connected_pop_2050'] = pipelines_NWC_edges_pop_update.groupby('node_id')['asset_connected_pop_2050'].sum().to_frame().reset_index()['asset_connected_pop_2050']
pipelines_NWC_edges_pop_update = pipelines_NWC_edges_pop_update_gb.merge(pipelines_NWC_edges_pop_update.drop(['asset_connected_pop_2020','asset_connected_pop_2030','asset_connected_pop_2040','asset_connected_pop_2050'],axis=1), on='node_id')
print(pipelines_NWC_edges_pop_update.groupby('WSZONEID')['asset_connected_pop_2020'].max().to_frame().reset_index()['asset_connected_pop_2020'].sum())
#pipelines_NWC_edges_pop_update = gpd.to_file('potable_facilities_NWC_pop_update_drop_dups.gpkg')
