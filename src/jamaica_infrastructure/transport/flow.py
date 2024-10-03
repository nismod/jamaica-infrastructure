"""Functions used in the provincial and national-scale network failure analysis
"""

from collections import defaultdict
from itertools import chain

import igraph as ig
import networkx as nx
import numpy as np
import pandas as pd
from tqdm import tqdm


def swap_min_max(x, min_col, max_col):
    """Swap columns if necessary"""
    if x[min_col] < 0 and x[max_col] < 0:
        if abs(x[min_col]) > abs(x[max_col]):
            return x[max_col], x[min_col]
        else:
            return x[min_col], x[max_col]
    else:
        if x[min_col] > x[max_col]:
            return x[max_col], x[min_col]
        else:
            return x[min_col], x[max_col]


def get_flow_on_edges(save_paths_df, edge_id_column, edge_path_column, flow_column):
    """Write results to Shapefiles

    Outputs ``gdf_edges`` - a shapefile with minimum and maximum tonnage flows of all
    commodities/industries for each edge of network.

    Parameters
    ---------
    save_paths_df
        Pandas DataFrame of OD flow paths and their tonnages
    industry_columns
        List of string names of all OD commodities/industries indentified
    min_max_exist
        List of string names of commodity/industry columns for which min-max tonnage column names already exist
    gdf_edges
        GeoDataFrame of network edge set
    save_csv
        Boolean condition to tell code to save created edge csv file
    save_shapes
        Boolean condition to tell code to save created edge shapefile
    shape_output_path
        Path where the output shapefile will be stored
    csv_output_path
        Path where the output csv file will be stored

    """
    edge_flows = defaultdict(float)
    for row in save_paths_df.itertuples():
        for item in getattr(row, edge_path_column):
            edge_flows[item] += getattr(row, flow_column)

    return pd.DataFrame(
        [(k, v) for k, v in edge_flows.items()], columns=[edge_id_column, flow_column]
    )


def get_flow_paths_indexes_of_edges(flow_dataframe, path_criteria):
    edge_path_index = defaultdict(list)
    for k, v in zip(
        chain.from_iterable(flow_dataframe[path_criteria].ravel()),
        flow_dataframe.index.repeat(flow_dataframe[path_criteria].str.len()).tolist(),
    ):
        edge_path_index[k].append(v)

    del flow_dataframe
    return edge_path_index


def get_flow_paths_numbers_of_edges(edge_path_index):
    edge_path_count = defaultdict(int)
    for k, v in edge_path_index.items():
        edge_path_count[k] = len(v)

    return edge_path_count


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


def igraph_scenario_edge_failures(
    network_df_in,
    edge_failure_set,
    flow_dataframe,
    edge_flow_path_indexes,
    path_criteria,
    cost_criteria,
    new_path=False,
):
    """Estimate network impacts of each failures
    When the tariff costs of each path are fixed by vehicle weight

    Parameters
    ---------
    network_df_in - Pandas DataFrame of network
    edge_failure_set - List of string edge ID's
    flow_dataframe - Pandas DataFrame of list of edge paths
    path_criteria - String name of column of edge paths in flow dataframe
    tons_criteria - String name of column of path tons in flow dataframe
    cost_criteria - String name of column of path costs in flow dataframe
    time_criteria - String name of column of path travel time in flow dataframe


    Returns
    -------
    edge_failure_dictionary : list[dict]
        With attributes
        edge_id - String name or list of failed edges
        origin - String node ID of Origin of disrupted OD flow
        destination - String node ID of Destination of disrupted OD flow
        no_access - Boolean 1 (no reroutng) or 0 (rerouting)
        new_cost - Float value of estimated cost of OD journey after disruption
        new_distance - Float value of estimated distance of OD journey after disruption
        new_path - List of string edge ID's of estimated new route of OD journey after disruption
        new_time - Float value of estimated time of OD journey after disruption
    """
    edge_fail_dictionary = []
    # network_df,edge_path_index = identify_all_failure_paths(network_df_in,edge_failure_set,flow_dataframe,path_criteria)

    edge_path_index = list(
        set(
            list(
                chain.from_iterable(
                    [
                        path_idx
                        for path_key, path_idx in edge_flow_path_indexes.items()
                        if path_key in edge_failure_set
                    ]
                )
            )
        )
    )

    if edge_path_index:
        select_flows = flow_dataframe[flow_dataframe.index.isin(edge_path_index)]
        del edge_path_index
        network_graph = ig.Graph.TupleList(
            network_df_in[~network_df_in["edge_id"].isin(edge_failure_set)].itertuples(
                index=False
            ),
            edge_attrs=list(network_df_in.columns)[2:],
        )

        first_edge_id = edge_failure_set[0]
        del edge_failure_set
        A = sorted(
            network_graph.clusters().subgraphs(),
            key=lambda l: len(l.es["edge_id"]),
            reverse=True,
        )
        access_flows = []
        edge_fail_dictionary = []
        for i in range(len(A)):
            network_graph = A[i]
            nodes_name = np.asarray([x["name"] for x in network_graph.vs])
            po_access = select_flows[
                (select_flows["origin_id"].isin(nodes_name))
                & (select_flows["destination_id"].isin(nodes_name))
            ]

            if len(po_access.index) > 0:
                po_access = po_access.set_index("origin_id")
                origins = list(set(po_access.index.values.tolist()))
                for o in range(len(origins)):
                    origin = origins[o]
                    destinations = po_access.loc[
                        [origin], "destination_id"
                    ].values.tolist()
                    # tons = po_access.loc[[origin], tons_criteria].values.tolist()
                    paths = network_graph.get_shortest_paths(
                        origin, destinations, weights=cost_criteria, output="epath"
                    )
                    if new_path is True:
                        for p in range(len(paths)):
                            new_gcost = 0
                            new_path = []
                            for n in paths[p]:
                                new_gcost += network_graph.es[n][cost_criteria]
                                new_path.append(network_graph.es[n]["edge_id"])
                            edge_fail_dictionary.append(
                                {
                                    "edge_id": first_edge_id,
                                    "origin_id": origin,
                                    "destination_id": destinations[p],
                                    "new_path": new_path,
                                    "new_cost": new_gcost,
                                    "no_access": 0,
                                }
                            )
                    else:
                        for p in range(len(paths)):
                            new_gcost = 0
                            for n in paths[p]:
                                new_gcost += network_graph.es[n][cost_criteria]
                            edge_fail_dictionary.append(
                                {
                                    "edge_id": first_edge_id,
                                    "origin_id": origin,
                                    "destination_id": destinations[p],
                                    "new_cost": new_gcost,
                                    "no_access": 0,
                                }
                            )
                    del destinations, paths
                del origins
                po_access = po_access.reset_index()
                po_access["access"] = 1
                access_flows.append(
                    po_access[["origin_id", "destination_id", "access"]]
                )
            del po_access

        del A

        if len(access_flows):
            access_flows = pd.concat(
                access_flows, axis=0, sort="False", ignore_index=True
            )
            select_flows = pd.merge(
                select_flows,
                access_flows,
                how="left",
                on=["origin_id", "destination_id"],
            ).fillna(0)
        else:
            select_flows["access"] = 0

        no_access = select_flows[select_flows["access"] == 0]
        if len(no_access.index) > 0:
            for value in no_access.itertuples():
                if new_path is True:
                    edge_fail_dictionary.append(
                        {
                            "edge_id": first_edge_id,
                            "origin_id": getattr(value, "origin_id"),
                            "destination_id": getattr(value, "destination_id"),
                            "new_path": [],
                            "new_cost": 0,
                            "no_access": 1,
                        }
                    )
                else:
                    edge_fail_dictionary.append(
                        {
                            "edge_id": first_edge_id,
                            "origin_id": getattr(value, "origin_id"),
                            "destination_id": getattr(value, "destination_id"),
                            "new_cost": 0,
                            "no_access": 1,
                        }
                    )

        del no_access, select_flows

    return edge_fail_dictionary


def rearrange_minmax_values(edge_failure_dataframe):
    """Write results to Shapefiles

    Parameters
    ---------
    edge_failure_dataframe : pandas.DataFrame
        with min-max columns

    Returns
    -------
    edge_failure_dataframe : pandas.DataFrame
        With columns where min < max
    """
    failure_columns = edge_failure_dataframe.columns.values.tolist()
    failure_columns = [f for f in failure_columns if f != ("edge_id", "no_access")]

    industry_columns = list(
        set([f.split("min_")[1] for f in failure_columns if "min" in f])
    )

    for ind in industry_columns:
        edge_failure_dataframe["swap"] = edge_failure_dataframe.apply(
            lambda x: swap_min_max(x, "min_{}".format(ind), "max_{}".format(ind)),
            axis=1,
        )
        edge_failure_dataframe[["min_{}".format(ind), "max_{}".format(ind)]] = (
            edge_failure_dataframe["swap"].apply(pd.Series)
        )
        edge_failure_dataframe.drop("swap", axis=1, inplace=True)

    return edge_failure_dataframe


def network_failure_assembly_shapefiles(
    edge_failure_dataframe, gdf_edges, save_edges=True, shape_output_path=""
):
    """Write results to Shapefiles


    Outputs gdf_edges - a Shapefile with results of edge failure dataframe

    Parameters
    ---------
    edge_failure_dataframe
        Pandas DataFrame of edge failure results
    gdf_edges
        GeoDataFrame of network edge set with edge ID's and geometry
    save_edges : bool
        Boolean condition to tell code to save created edge shapefile
    shape_output_path : str
        Path where the output shapefile will be stored

    """
    failure_columns = edge_failure_dataframe.columns.values.tolist()
    failure_columns = [f for f in failure_columns if f != "edge_id"]

    for fc in failure_columns:
        gdf_edges[fc] = 0

    for iter_, row in edge_failure_dataframe.iterrows():
        # print (row[1:])
        gdf_edges.loc[gdf_edges["edge_id"] == row["edge_id"], failure_columns] = row[
            failure_columns
        ].values

    industry_columns = list(
        set([f.split("min_")[1] for f in failure_columns if "min" in f])
    )

    for ind in industry_columns:
        gdf_edges["swap"] = gdf_edges.apply(
            lambda x: swap_min_max(x, "min_{}".format(ind), "max_{}".format(ind)),
            axis=1,
        )
        gdf_edges[["min_{}".format(ind), "max_{}".format(ind)]] = gdf_edges[
            "swap"
        ].apply(pd.Series)
        gdf_edges.drop("swap", axis=1, inplace=True)

    if save_edges == True:
        gdf_edges.to_file(shape_output_path)

    del gdf_edges, edge_failure_dataframe


def edge_failure_sampling(failure_scenarios, edge_column):
    """Criteria for selecting failure samples

    Parameters
    ---------
    failure_scenarios - Pandas DataFrame of failure scenarios
    edge_column - String name of column to select failed edge ID's

    Returns
    -------
    edge_failure_samples - List of lists of failed edge sets
    """
    edge_failure_samples = list(set(failure_scenarios[edge_column].values.tolist()))

    return edge_failure_samples


def merge_failure_results(
    flow_df_select, failure_df, id_col, tons_col, dist_col, time_col, cost_col
):
    """Merge failure results with flow results

    Parameters
    ---------
    flow_df_select : pandas.DataFrame
        edge flow values
    failure_df : pandas.DataFrame
        edge failure values
    tons_col : str
        name of column of tonnages in flow dataframe
    dist_col : str
        name of column of distance in flow dataframe
    time_col : str
        name of column of time in flow dataframe
    cost_col : str
        name of column of cost in flow dataframe
    vehicle_col : str
        name of column of vehicle counts in flow dataframe
    changing_tonnages : bool

    Returns
    -------
    flow_df_select : pandas.DataFrame
        Of edge flow and failure values merged
    """
    flow_df_select = pd.merge(
        flow_df_select, failure_df, on=["origin_id", "destination_id"], how="left"
    ).fillna(0)
    flow_df_select = flow_df_select[
        (flow_df_select[tons_col] > 0) & (flow_df_select[id_col] != 0)
    ]

    flow_df_select["dist_diff"] = (1 - flow_df_select["no_access"]) * (
        flow_df_select["new_distance"] - flow_df_select[dist_col]
    )
    flow_df_select["time_diff"] = (1 - flow_df_select["no_access"]) * (
        flow_df_select["new_time"] - flow_df_select[time_col]
    )
    flow_df_select["tr_loss"] = (1 - flow_df_select["no_access"]) * (
        flow_df_select["new_cost"] - flow_df_select[cost_col]
    )

    return flow_df_select
