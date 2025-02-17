"""Functions for preprocessing data
"""

import sys
import os
import json

import pandas as pd
import geopandas as gpd
from scipy.spatial import Voronoi
from shapely.geometry import Polygon, shape, LineString

# workaround for geopandas >0.9 until snkit #37 and geopandas #1977 are fixed
gpd._compat.USE_PYGEOS = False
import fiona
import numpy as np
from scipy.spatial import cKDTree
import snkit
from tqdm import tqdm

tqdm.pandas()
epsg_jamaica = 3448


def load_config():
    """Read config.json"""
    config_path = os.path.join(os.path.dirname(__file__), "..", "..", "config.json")
    with open(config_path, "r") as config_fh:
        config = json.load(config_fh)
    return config


def geopandas_read_file_type(file_path, file_layer, file_database=None):
    if file_database is not None:
        return gpd.read_file(os.path.join(file_path, file_database), layer=file_layer)
    else:
        return gpd.read_file(os.path.join(file_path, file_layer))


def nearest(geom, gdf):
    """Find the element of a GeoDataFrame nearest a shapely geometry"""
    matches_idx = gdf.sindex.nearest(geom.bounds)
    nearest_geom = min(
        [gdf.iloc[match_idx] for match_idx in matches_idx],
        key=lambda match: geom.distance(match.geometry),
    )
    return nearest_geom


def get_nearest_node(x, sindex_input_nodes, input_nodes, id_column):
    """Get nearest node in a dataframe

    Parameters
    ----------
    x
        row of dataframe
    sindex_nodes
        spatial index of dataframe of nodes in the network
    nodes
        dataframe of nodes in the network
    id_column
        name of column of id of closest node

    Returns
    -------
    Nearest node to geometry of row
    """
    return input_nodes.loc[list(sindex_input_nodes.nearest(x.bounds[:2]))][
        id_column
    ].values[0]


# def add_topology(network, id_col="id"):
#     """Add or replace from_id, to_id to edges"""
#     from_ids = []
#     to_ids = []

#     for edge in tqdm(
#         network.edges.itertuples(), desc="topology", total=len(network.edges)
#     ):
#         start, end = snkit.network.line_endpoints(edge.geometry)

#         start_node = nearest(start, network.nodes)
#         from_ids.append(start_node[id_col])

#         end_node = nearest(end, network.nodes)
#         to_ids.append(end_node[id_col])

#     edges = network.edges.copy()
#     edges["from_id"] = from_ids
#     edges["to_id"] = to_ids

#     return snkit.network.Network(nodes=network.nodes, edges=edges)


def create_network_from_nodes_and_edges(
    nodes, edges, node_edge_prefix, snap_distance=None, by=None
):
    edges.columns = map(str.lower, edges.columns)
    if "id" in edges.columns.values.tolist():
        edges.rename(columns={"id": "e_id"}, inplace=True)

    # Deal with empty edges (drop)
    empty_idx = edges.geometry.apply(lambda e: e is None or e.is_empty)
    if empty_idx.sum():
        empty_edges = edges[empty_idx]
        print(f"Found {len(empty_edges)} empty edges.")
        print(empty_edges)
        edges = edges[~empty_idx].copy()

    network = snkit.Network(nodes, edges)
    print("* Done with network creation")

    network = snkit.network.split_multilinestrings(network)
    print("* Done with splitting multilines")

    if nodes is not None:
        if snap_distance is not None:
            network = snkit.network.link_nodes_to_edges_within(
                network, snap_distance, tolerance=1e-09
            )
            print("* Done with joining nodes to edges")
        else:
            network = snkit.network.snap_nodes(network)
            print("* Done with snapping nodes to edges")
        # network.nodes = snkit.network.drop_duplicate_geometries(network.nodes)
        # print ('* Done with dropping same geometries')

        # network = snkit.network.split_edges_at_nodes(network,tolerance=9e-10)
        # print ('* Done with splitting edges at nodes')

    network = snkit.network.add_endpoints(network)
    print("* Done with adding endpoints")

    network.nodes = snkit.network.drop_duplicate_geometries(network.nodes)
    print("* Done with dropping same geometries")

    network = snkit.network.split_edges_at_nodes(network, tolerance=1e-3)
    print("* Done with splitting edges at nodes")

    network = snkit.network.add_ids(
        network, edge_prefix=f"{node_edge_prefix}e", node_prefix=f"{node_edge_prefix}n"
    )
    network = snkit.network.add_topology(network, id_col="id")
    print("* Done with network topology")

    if by is not None:
        network = snkit.network.merge_edges(network, by=by)
        print("* Done with merging network")

    network.edges.rename(
        columns={"from_id": "from_node", "to_id": "to_node", "id": "edge_id"},
        inplace=True,
    )
    network.nodes.rename(columns={"id": "node_id"}, inplace=True)

    return network


def add_network_topology(nodes, edges, node_edge_prefix):
    network = snkit.Network(nodes=nodes, edges=edges)
    network = snkit.network.add_ids(
        network, edge_prefix=f"{node_edge_prefix}e", node_prefix=f"{node_edge_prefix}n"
    )
    network = snkit.network.add_topology(network, id_col="id")
    print("* Done with network topology")

    network = snkit.network.add_component_ids(network)
    print("* Done with add network components")

    network.edges.rename(
        columns={"from_id": "from_node", "to_id": "to_node", "id": "edge_id"},
        inplace=True,
    )
    network.nodes.rename(columns={"id": "node_id"}, inplace=True)

    return network


def voronoi_finite_polygons_2d(vor, radius=None):
    """Reconstruct infinite voronoi regions in a 2D diagram to finite regions.

    Source: https://stackoverflow.com/questions/36063533/clipping-a-voronoi-diagram-python

    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'

    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end
    """
    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max() * 2

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1]  # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:, 1] - c[1], vs[:, 0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)


def assign_value_in_area_proportions(poly_1_gpd, poly_2_gpd, poly_attribute):
    poly_1_sindex = poly_1_gpd.sindex
    for p_2_index, polys_2 in poly_2_gpd.iterrows():
        poly2_attr = 0
        intersected_polys = poly_1_gpd.iloc[
            list(poly_1_sindex.intersection(polys_2.geometry.bounds))
        ]
        for p_1_index, polys_1 in intersected_polys.iterrows():
            if (
                (polys_2["geometry"].intersects(polys_1["geometry"]) is True)
                and (polys_1.geometry.is_valid is True)
                and (polys_2.geometry.is_valid is True)
            ):
                poly2_attr += (
                    polys_1[poly_attribute]
                    * polys_2["geometry"].intersection(polys_1["geometry"]).area
                    / polys_1["geometry"].area
                )

        poly_2_gpd.loc[p_2_index, poly_attribute] = poly2_attr

    return poly_2_gpd


def extract_nodes_within_gdf(x, input_nodes, column_name):
    a = input_nodes.loc[list(input_nodes.geometry.within(x.geometry))]
    # if len(a.index) > 1: # To check if there are multiple intersections
    #     print (x)
    if len(a.index) > 0:
        return a[column_name].values[0]
    else:
        return ""


def create_voronoi_layer(nodes_dataframe, node_id_column, epsg=4326, **kwargs):
    """Assign weights to nodes based on their nearest populations

        - By finding the population that intersect with the Voronoi extents of nodes

    Parameters
        - nodes_dataframe - Geodataframe of the nodes
        - population_dataframe - Geodataframe of the population
        - nodes_id_column - String name of node ID column
        - population_value_column - String name of column containing population values

    Outputs
        - nodes - Geopandas dataframe of nodes with new column called population
    """

    # load provinces and get geometry of the right population_dataframe

    # create Voronoi polygons for the nodes
    xy_list = []
    for iter_, values in nodes_dataframe.iterrows():
        xy = list(values.geometry.coords)
        xy_list += [list(xy[0])]

    vor = Voronoi(np.array(xy_list))
    regions, vertices = voronoi_finite_polygons_2d(vor)
    min_x = vor.min_bound[0] - 0.1
    max_x = vor.max_bound[0] + 0.1
    min_y = vor.min_bound[1] - 0.1
    max_y = vor.max_bound[1] + 0.1

    mins = np.tile((min_x, min_y), (vertices.shape[0], 1))
    bounded_vertices = np.max((vertices, mins), axis=0)
    maxs = np.tile((max_x, max_y), (vertices.shape[0], 1))
    bounded_vertices = np.min((bounded_vertices, maxs), axis=0)

    box = Polygon([[min_x, min_y], [min_x, max_y], [max_x, max_y], [max_x, min_y]])

    poly_list = []
    for region in regions:
        polygon = vertices[region]
        print(region)
        print(polygon)
        # Clipping polygon
        poly = Polygon(polygon)
        poly = poly.intersection(box)
        poly_list.append(poly)

    poly_index = list(np.arange(0, len(poly_list), 1))
    poly_df = pd.DataFrame(
        list(zip(poly_index, poly_list)), columns=["gid", "geometry"]
    )
    gdf_voronoi = gpd.GeoDataFrame(poly_df, geometry="geometry", crs=f"epsg:{epsg}")
    gdf_voronoi["areas"] = gdf_voronoi.progress_apply(lambda x: x.geometry.area, axis=1)
    gdf_voronoi[node_id_column] = gdf_voronoi.progress_apply(
        lambda x: extract_nodes_within_gdf(x, nodes_dataframe, node_id_column), axis=1
    )
    if not kwargs.get("save", False):
        pass
    else:
        gdf_voronoi.to_file(kwargs.get("voronoi_path", "voronoi-output.shp"))

    return gdf_voronoi


def assign_node_weights_by_population_proximity(
    nodes_dataframe,
    population_dataframe,
    node_id_column,
    population_value_column,
    epsg=4326,
    **kwargs,
):
    """Assign weights to nodes based on their nearest populations

        - By finding the population that intersect with the Voronoi extents of nodes

    Parameters
        - nodes_dataframe - Geodataframe of the nodes
        - population_dataframe - Geodataframe of the population
        - nodes_id_column - String name of node ID column
        - population_value_column - String name of column containing population values

    Outputs
        - nodes - Geopandas dataframe of nodes with new column called population
    """

    # load provinces and get geometry of the right population_dataframe
    sindex_population_dataframe = population_dataframe.sindex

    # create Voronoi polygons for the nodes
    xy_list = []
    for iter_, values in nodes_dataframe.iterrows():
        xy = list(values.geometry.coords)
        xy_list += [list(xy[0])]

    vor = Voronoi(np.array(xy_list))
    regions, vertices = voronoi_finite_polygons_2d(vor)
    min_x = vor.min_bound[0] - 0.1
    max_x = vor.max_bound[0] + 0.1
    min_y = vor.min_bound[1] - 0.1
    max_y = vor.max_bound[1] + 0.1

    mins = np.tile((min_x, min_y), (vertices.shape[0], 1))
    bounded_vertices = np.max((vertices, mins), axis=0)
    maxs = np.tile((max_x, max_y), (vertices.shape[0], 1))
    bounded_vertices = np.min((bounded_vertices, maxs), axis=0)

    box = Polygon([[min_x, min_y], [min_x, max_y], [max_x, max_y], [max_x, min_y]])

    poly_list = []
    for region in regions:
        polygon = vertices[region]
        # Clipping polygon
        poly = Polygon(polygon)
        poly = poly.intersection(box)
        poly_list.append(poly)

    poly_index = list(np.arange(0, len(poly_list), 1))
    poly_df = pd.DataFrame(
        list(zip(poly_index, poly_list)), columns=["gid", "geometry"]
    )
    gdf_voronoi = gpd.GeoDataFrame(poly_df, geometry="geometry", crs=f"epsg:{epsg}")
    gdf_voronoi["areas"] = gdf_voronoi.progress_apply(lambda x: x.geometry.area, axis=1)
    gdf_voronoi[node_id_column] = gdf_voronoi.progress_apply(
        lambda x: extract_nodes_within_gdf(x, nodes_dataframe, node_id_column), axis=1
    )
    if not kwargs.get("save", False):
        pass
    else:
        gdf_voronoi.to_file(kwargs.get("voronoi_path", "voronoi-output.shp"))

    gdf_voronoi[population_value_column] = 0
    gdf_voronoi = assign_value_in_area_proportions(
        population_dataframe, gdf_voronoi, population_value_column
    )
    gdf_voronoi = gdf_voronoi[~(gdf_voronoi[node_id_column] == "")]

    gdf_pops = gdf_voronoi[[node_id_column, population_value_column]]
    del gdf_voronoi, poly_list, poly_df

    nodes_dataframe = pd.merge(
        nodes_dataframe, gdf_pops, how="left", on=[node_id_column]
    )
    del gdf_pops, population_dataframe

    return nodes_dataframe


def spatial_scenario_selection(
    dataframe_1,
    dataframe_2,
    dataframe_1_columns,
    dataframe_2_columns,
):
    """Intersect Polygons to collect attributes

    Parameters
        - dataframe_1 - First polygon dataframe
        - dataframe_2 - Second polygon dataframe
        - dataframe_1_columns - First polygon dataframe columns to collect
        - dataframe_2_columns - Second polygon dataframe columns to collect

    Outputs
        data_dictionary - Dictionary of intersection attributes:
    """

    intersection_dictionary = []

    # create spatial index
    dataframe_1_sindex = dataframe_1.sindex
    total_values = len(dataframe_2.index)
    for values in dataframe_2.itertuples():
        intersected_polys = dataframe_1.iloc[
            list(dataframe_1_sindex.intersection(values.geometry.bounds))
        ]
        for intersected_values in intersected_polys.itertuples():
            if (
                (intersected_values.geometry.intersects(values.geometry) is True)
                and (values.geometry.is_valid is True)
                and (intersected_values.geometry.is_valid is True)
            ):
                dataframe_1_dictionary = dict(
                    [(v, getattr(intersected_values, v)) for v in dataframe_1_columns]
                )
                dataframe_2_dictionary = dict(
                    [(v, getattr(values, v)) for v in dataframe_2_columns]
                )
                geometry_dictionary = {
                    "geometry": values.geometry.intersection(
                        intersected_values.geometry
                    )
                }

                intersection_dictionary.append(
                    {
                        **dataframe_1_dictionary,
                        **dataframe_2_dictionary,
                        **geometry_dictionary,
                    }
                )
        print(f"* Done with Index {values.Index} out of {total_values}")
    return intersection_dictionary


def split_multigeometry(dataframe, split_geometry_type="GeometryCollection"):
    """Create multiple geometries from any MultiGeomtery and GeometryCollection

    Ensures that edge geometries are all Points,LineStrings,Polygons, duplicates attributes over any
    created multi-geomteries.
    """
    simple_geom_attrs = []
    simple_geom_geoms = []
    for v in tqdm(
        dataframe.itertuples(index=False), desc="split_multi", total=len(dataframe)
    ):
        if v.geometry.geom_type == split_geometry_type:
            geom_parts = list(v.geometry)
        else:
            geom_parts = [v.geometry]

        for part in geom_parts:
            simple_geom_geoms.append(part)

        attrs = gpd.GeoDataFrame([v] * len(geom_parts))
        simple_geom_attrs.append(attrs)

    simple_geom_geoms = gpd.GeoDataFrame(simple_geom_geoms, columns=["geometry"])
    dataframe = (
        pd.concat(simple_geom_attrs, axis=0)
        .reset_index(drop=True)
        .drop("geometry", axis=1)
    )
    dataframe = pd.concat([dataframe, simple_geom_geoms], axis=1)

    return dataframe


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
    nearest_pts["geometry"] = nearest_pts.progress_apply(
        lambda x: LineString([x.from_geometry, x.to_geometry]), axis=1
    )
    nearest_pts.drop(["from_geometry", "to_geometry"], axis=1, inplace=True)
    nearest_pts["from_mode"] = from_mode
    nearest_pts["to_mode"] = to_mode
    nearest_pts = gpd.GeoDataFrame(
        nearest_pts, geometry="geometry", crs=f"EPSG:{epsg_jamaica}"
    )

    return nearest_pts
