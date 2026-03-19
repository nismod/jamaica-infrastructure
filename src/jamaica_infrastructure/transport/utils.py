from collections import defaultdict
import logging

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from shapely.geometry import LineString
from tqdm import tqdm

from jamaica_infrastructure.geo import LOCAL_PROJ_CRS_EPSG


# workaround for geopandas >0.9 until snkit #37 and geopandas #1977 are fixed
gpd._compat.USE_PYGEOS = False

tqdm.pandas()


def get_flow_on_edges(save_paths_df, edge_id_column, edge_path_column, flow_column):
    edge_flows = defaultdict(float)
    for row in save_paths_df.itertuples():
        for item in getattr(row, edge_path_column):
            edge_flows[item] += getattr(row, flow_column)

    return pd.DataFrame(
        [(k, v) for k, v in edge_flows.items()], columns=[edge_id_column, flow_column]
    )


def nearest(geom, gdf):
    """Find the element of a GeoDataFrame nearest a shapely geometry"""
    matches_idx = gdf.sindex.nearest(geom.bounds)
    nearest_geom = min(
        [gdf.iloc[match_idx] for match_idx in matches_idx],
        key=lambda match: geom.distance(match.geometry),
    )
    return nearest_geom


def network_od_path_estimations(graph, source, target, cost_criteria):
    """Estimate the paths, distances, times, and costs for given OD pair

    Parameters
    ---------
    graph
        igraph network structure
    source
        String/Float/Integer name of Origin node ID
    source
        String/Float/Integer name of Destination node ID
    tonnage : float
        value of tonnage
    vehicle_weight : float
        unit weight of vehicle
    cost_criteria : str
        name of generalised cost criteria to be used: min_gcost or max_gcost
    time_criteria : str
        name of time criteria to be used: min_time or max_time
    fixed_cost : bool

    Returns
    -------
    edge_path_list : list[list]
        nested lists of Strings/Floats/Integers of edge ID's in routes
    path_dist_list : list[float]
        estimated distances of routes
    path_time_list : list[float]
        estimated times of routes
    path_gcost_list : list[float]
        estimated generalised costs of routes

    """
    paths = graph.get_shortest_paths(
        source, target, weights=cost_criteria, output="epath"
    )

    edge_path_list = []
    path_gcost_list = []
    # for p in range(len(paths)):
    for path in paths:
        edge_path = []
        path_gcost = 0
        if path:
            for n in path:
                edge_path.append(graph.es[n]["edge_id"])
                path_gcost += graph.es[n][cost_criteria]

        edge_path_list.append(edge_path)
        path_gcost_list.append(path_gcost)

    return edge_path_list, path_gcost_list


def network_od_paths_assembly(points_dataframe, graph, cost_criteria, tonnage_column):
    """Assemble estimates of OD paths, distances, times, costs and tonnages on networks

    Parameters
    ----------
    points_dataframe : pandas.DataFrame
        OD nodes and their tonnages
    graph
        igraph network structure
    region_name : str
        name of Province
    excel_writer
        Name of the excel writer to save Pandas dataframe to Excel file

    Returns
    -------
    save_paths_df : pandas.DataFrame
        - origin - String node ID of Origin
        - destination - String node ID of Destination
        - min_edge_path - List of string of edge ID's for paths with minimum generalised cost flows
        - max_edge_path - List of string of edge ID's for paths with maximum generalised cost flows
        - min_netrev - Float values of estimated netrevenue for paths with minimum generalised cost flows
        - max_netrev - Float values of estimated netrevenue for paths with maximum generalised cost flows
        - min_croptons - Float values of estimated crop tons for paths with minimum generalised cost flows
        - max_croptons - Float values of estimated crop tons for paths with maximum generalised cost flows
        - min_distance - Float values of estimated distance for paths with minimum generalised cost flows
        - max_distance - Float values of estimated distance for paths with maximum generalised cost flows
        - min_time - Float values of estimated time for paths with minimum generalised cost flows
        - max_time - Float values of estimated time for paths with maximum generalised cost flows
        - min_gcost - Float values of estimated generalised cost for paths with minimum generalised cost flows
        - max_gcost - Float values of estimated generalised cost for paths with maximum generalised cost flows

    """
    save_paths = []
    points_dataframe = points_dataframe.set_index("origin_id")
    origins = list(set(points_dataframe.index.values.tolist()))
    for origin in origins:
        try:
            destinations = points_dataframe.loc[
                [origin], "destination_id"
            ].values.tolist()

            get_path, get_gcost = network_od_path_estimations(
                graph, origin, destinations, cost_criteria
            )

            # tons = points_dataframe.loc[[origin], tonnage_column].values
            save_paths += list(
                zip([origin] * len(destinations), destinations, get_path, get_gcost)
            )

        except:
            logging.info(f"* no path between {origin}-{destinations}")

    cols = ["origin_id", "destination_id", "edge_path", "gcost"]
    save_paths_df = pd.DataFrame(save_paths, columns=cols)
    points_dataframe = points_dataframe.reset_index()
    save_paths_df = pd.merge(
        save_paths_df, points_dataframe, how="left", on=["origin_id", "destination_id"]
    ).fillna(0)

    save_paths_df = save_paths_df[
        (save_paths_df[tonnage_column] > 0) & (save_paths_df["origin_id"] != 0)
    ]

    return save_paths_df


def ckdnearest(gdA, gdB):
    """Taken from https://gis.stackexchange.com/questions/222315/finding-nearest-point-in-other-geodataframe-using-geopandas"""
    nA = np.array(list(gdA.geometry.apply(lambda x: (x.x, x.y))))
    nB = np.array(list(gdB.geometry.apply(lambda x: (x.x, x.y))))
    btree = cKDTree(nB)
    dist, idx = btree.query(nA, k=1)
    gdB_nearest = gdB.iloc[idx].drop(columns="geometry").reset_index(drop=True)
    gdf = pd.concat(
        [gdA.reset_index(drop=True), gdB_nearest, pd.Series(dist, name="dist")], axis=1
    )

    return gdf


def polygon_to_points(gdf):
    gdf = gdf.to_crs(epsg=LOCAL_PROJ_CRS_EPSG)
    gdf["geometry"] = gdf.apply(lambda x: x.geometry.centroid, axis=1)

    return gdf


def map_nearest_locations_and_create_lines(
    from_gdf, to_gdf, from_gdf_id, to_gdf_id, from_mode, to_mode
):
    from_gdf.rename(columns={from_gdf_id: "from_node"}, inplace=True)
    to_gdf.rename(columns={to_gdf_id: "to_node"}, inplace=True)
    nearest_pts = ckdnearest(
        from_gdf[["from_node", "geometry"]], to_gdf[["to_node", "geometry"]]
    )
    nearest_pts.rename(
        columns={"geometry": "from_geometry", "dist": "length_m"}, inplace=True
    )
    nearest_pts = pd.merge(
        nearest_pts, to_gdf[["to_node", "geometry"]], how="left", on=["to_node"]
    )
    nearest_pts.rename(columns={"geometry": "to_geometry"}, inplace=True)
    nearest_pts["geometry"] = nearest_pts.apply(
        lambda x: LineString([x.from_geometry, x.to_geometry]), axis=1
    )
    nearest_pts.drop(["from_geometry", "to_geometry"], axis=1, inplace=True)
    nearest_pts["from_mode"] = from_mode
    nearest_pts["to_mode"] = to_mode
    nearest_pts = gpd.GeoDataFrame(
        nearest_pts, geometry="geometry", crs=f"EPSG:{LOCAL_PROJ_CRS_EPSG}"
    )

    return nearest_pts
