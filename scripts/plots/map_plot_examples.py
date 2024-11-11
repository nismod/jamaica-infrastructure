# shapefiles
island = gpd.read_file(
    "/soge-home/projects/mistral/jamaica-ccri/incoming_data/water/potable/jamaica-latest-free.shp/jam_admbnda_adm0_regrid.shp",
    crs="EPSG:3448",
)
pipeline = gpd.read_file(
    "/soge-home/projects/mistral/jamaica-ccri/incoming_data/water/potable/raw/pipelines_network_NWC.shp",
    crs="EPSG:3448",
)
river_network = gpd.read_file(
    "/soge-home/projects/mistral/jamaica-ccri/incoming_data/water/natural/river_network.shp",
    crs="EPSG:3448",
)
major_rivers = gpd.read_file(
    "/soge-home/projects/mistral/jamaica-ccri/incoming_data/water/natural/major_rivers.shp",
    crs="EPSG:3448",
)
hydrobasins = gpd.read_file(
    "/soge-home/projects/mistral/jamaica-ccri/incoming_data/water/natural/hydrobasins.shp",
    crs="EPSG:3448",
)

licensed_sources = gpd.read_file(
    "/soge-home/projects/mistral/jamaica-ccri/drought/economic/licensed_abstractions.shp",
    crs="EPSG:3448",
)
licensed_sources["lon"] = licensed_sources["geometry"].centroid.x
licensed_sources["lat"] = licensed_sources["geometry"].centroid.y
# licensed_sources = pd.read_csv('licensed sources submanagement catchment 2.csv')[['hybnnum', 'OBJECTID_2', 'Expiration', 'Source_Typ', 'Purpose_Us', 'VBF_Class', 'Lic_Stat','VolGrantm3']]
licensed_sources = licensed_sources[
    licensed_sources["Expiration"].isin(["Not Expired", "#VALUE!", ""])
]
licensed_sources["Source"] = np.where(
    licensed_sources["Source_Typ"] == "Well", "GW", "SW"
)
licensed_sources["Purpose_Us"] = licensed_sources["Purpose_Us"].astype(str)

licensed_sources["VBF_Class"] = np.where(
    licensed_sources["Purpose_Us"].str.contains("Public"),
    "Potable",
    np.where(
        licensed_sources["Purpose_Us"].str.contains("Potable"),
        "Potable",
        np.where(
            licensed_sources["Purpose_Us"].str.contains("Irrigation"),
            "Irrigation",
            np.where(
                licensed_sources["VBF_Class"].str.contains("Agricultural"),
                "Agriculture",
                np.where(
                    licensed_sources["Purpose_Us"].str.contains("Indus"),
                    "Industrial",
                    licensed_sources["VBF_Class"],
                ),
            ),
        ),
    ),
)

licensed_sources_sw = licensed_sources[licensed_sources["Source"] == "SW"]
licensed_sources_gw = licensed_sources[licensed_sources["Source"] == "GW"]
print(licensed_sources_sw)
# print(oliv)
fig, ax = plt.subplots(1, 1, figsize=(15, 7.5))
island.plot(ax=ax, color="none", edgecolor="black", linewidth=0.5)
groups = licensed_sources_sw.groupby("VBF_Class")
dict = {
    "Agriculture": "green",
    "Domestic": "blue",
    "Hydropower": "navy",
    "Industrial": "darkgoldenrod",
    "Irrigation": "lime",
    "Potable": "cyan",
    "Recreation": "purple",
}
for name, group in groups:
    c = dict[name]
    ax.scatter(
        group.lon,
        group.lat,
        marker="o",
        c=c,
        label=name + "SW",
        s=group.VolGrantm3 / 250,
        edgecolor="k",
    )  # , edgecolor='k'

groups = licensed_sources_gw.groupby("VBF_Class")
dict = {
    "Agriculture": "green",
    "Domestic": "blue",
    "Hydropower": "navy",
    "Industrial": "darkgoldenrod",
    "Irrigation": "lime",
    "Potable": "cyan",
    "Recreation": "purple",
}
for name, group in groups:
    c = dict[name]
    ax.scatter(
        group.lon,
        group.lat,
        marker="d",
        c=c,
        label=name + "GW",
        s=group.VolGrantm3 / 250,
        edgecolor="k",
    )  #

msizes = [
    round(licensed_sources["VolGrantm3"].quantile(0.25) / 250, 0),
    round(licensed_sources["VolGrantm3"].quantile(0.5) / 250, 0),
    round(licensed_sources["VolGrantm3"].quantile(0.75) / 250, 0),
    round(licensed_sources["VolGrantm3"].quantile(0.95) / 250, 0),
]
markers = []
for size in msizes:
    markers.append(
        ax.scatter(
            [],
            [],
            c="k",
            s=size,
            label=str(size * 500) + " Annual volume of abstraction (m3)",
        )
    )
lgnd = ax.legend()
lgnd.legendHandles[0]._sizes = [20]
lgnd.legendHandles[1]._sizes = [20]
lgnd.legendHandles[2]._sizes = [20]
lgnd.legendHandles[3]._sizes = [20]
lgnd.legendHandles[4]._sizes = [20]
lgnd.legendHandles[5]._sizes = [20]
lgnd.legendHandles[6]._sizes = [20]
lgnd.legendHandles[7]._sizes = [20]
lgnd.legendHandles[8]._sizes = [20]
lgnd.legendHandles[9]._sizes = [20]
lgnd.legendHandles[10]._sizes = [20]
lgnd.legendHandles[11]._sizes = [20]
lgnd.legendHandles[12]._sizes = [20]
lgnd.legendHandles[13]._sizes = [20]
ax.set_ylim(600000, 720000)
ax.set_xlim(600000, 850000)
plt.savefig("licensed_abstractions_4.png", dpi=300)

# # maps
# fig, ax = plt.subplots(2, 1, figsize=(20,30))
# island.plot(ax=ax[0], color='none', edgecolor='black', linewidth=0.5)

# hydrobasins_st = hydrobasins.dropna(subset=['hybnnum'])
# hydrobasins_st['lon'] = hydrobasins_st['geometry'].centroid.x
# hydrobasins_st['lat'] = hydrobasins_st['geometry'].centroid.y
# # for it,row in hydrobasins_st.iterrows():
# # ax[0].annotate(str(row['hybnnum']), (row['lon'], row['lat']))
# ax[0].scatter(filtered_streamflow_gauges_coords['POINT_X'],filtered_streamflow_gauges_coords['POINT_Y'], c='darkgrey', edgecolor='k', label='streamflow gauges (> 30 years data)')
# river_network.plot(ax=ax[0], color='lightblue', edgecolor='mediumblue', linewidth=0.25, label='rivers')
# major_rivers.plot(ax=ax[0], color='lightblue', edgecolor='mediumblue', linewidth=0.5)
# hydrobasins.plot(ax=ax[0], color='none', edgecolor='black', linewidth=0.5, label='catchments')
# ax[0].legend(fontsize=15)
# ax[0].set_ylim(600000,720000)
# ax[0].set_xlim(600000,850000)
# # ax[0].set_yticklabels(fontsize=20)
# # ax[0].set_xticklabels(fontsize=20)
# ax[0].tick_params(axis='both', which='major', labelsize=20)
# # ax[0].set_title('River network and streamflow gauges with >30 years of data')

# groups = potable_facilities_NWC.groupby('asset_type_cost_data')
# dict = {'treatment plant':'orangered', 'pumping unit':'grey', 'storage tank':'navy', 'well':'darkgoldenrod', 'intake':'cyan', }
# for name, group in groups:
# c = dict[name]
# ax[1].scatter(group.lon, group.lat, marker='o', c=c, label=name) #, edgecolor='k'
# pipeline.plot(ax=ax[1], color='black', edgecolor='black', linewidth=0.5, label='pipelines')
# island.plot(ax=ax[1], color='none', edgecolor='black', linewidth=0.5)
# ax[1].set_ylim(600000,720000)
# ax[1].set_xlim(600000,850000)
# # ax[1].set_yticklabels(fontsize=20)
# # ax[1].set_xticklabels(fontsize=20)
# ax[1].tick_params(axis='both', which='major', labelsize=20)
# ax[1].legend(fontsize=15)
# # ax[1].set_title('Potable water supply assets') #-76.42, -76.31 ##
# plt.savefig('network.png', dpi = 300)


# fig, ax = plt.subplots(1, 1, figsize=(20,10))

# buildings = gpd.read_file('jamaica-latest-free.shp/gis_osm_buildings_a_free_1_regrid.shp', crs="EPSG:3448")
# buildings.plot(ax=ax, color='lightgrey', edgecolor='lightgrey', linewidth=0.5, label='buildings')
# island.plot(ax=ax, color='none', edgecolor='black', linewidth=0.5)
# pipeline_convex.plot(ax=ax, color='none', edgecolor='red', linestyle='--', linewidth=0.75, label='system service areas')
# ax.scatter(potable_facilities_NWC.lon, potable_facilities_NWC.lat, marker='o', c='k', alpha=0.5, s=potable_facilities_NWC.asset_pop_new/500, edgecolor='k')
# msizes = [round(potable_facilities_NWC['asset_pop_new'].quantile(0.25)/500,0),
# round(potable_facilities_NWC['asset_pop_new'].quantile(0.5)/500,0),
# round(potable_facilities_NWC['asset_pop_new'].quantile(0.75)/500,0),
# round(potable_facilities_NWC['asset_pop_new'].quantile(0.95)/500,0)]
# markers = []
# for size in msizes:
# markers.append(ax.scatter([],[], c='k', s=size, label=str(size*500)+' people served per asset'))
# scen_85_patch = mpatches.Patch(color='grey', label='buildings')
# markers = markers.append(scen_85_patch)
# ax.legend(handles=markers,fontsize=15)
# ax.set_ylim(600000,720000) #16.95, 16.86
# ax.set_xlim(600000,850000) #-76.42, -76.31
# # ax.set_yticklabels(fontsize=20)
# # ax.set_xticklabels(fontsize=20)
# ax.tick_params(axis='both', which='major', labelsize=20)
# # ax.set_title('Asset criticality and system service areas')
# plt.savefig('network_2.png', dpi = 300)
