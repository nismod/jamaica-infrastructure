from collections import defaultdict
from itertools import chain
import json
import logging
import os

import geopandas as gpd
import igraph as ig
import numpy as np
import pandas as pd


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


def get_flow_paths_indexes_of_edges(flow_dataframe, path_criteria):
    edge_path_index = defaultdict(list)
    for k, v in zip(
        chain.from_iterable(flow_dataframe[path_criteria].to_numpy()),
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


def igraph_scenario_edge_failures_premade_network(
    network_graph: ig.Graph,
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
    network_graph - Premade iGraph network to disrupt
    edge_failure_set - List of string edge ID's
    flow_dataframe - Pandas DataFrame of list of edge paths
    path_criteria - String name of column of edge paths in flow dataframe
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

    for edge_id in edge_failure_set:
        try:
            network_graph.es.find(edge_id=edge_id).delete()
        except ValueError as error:
            if "no such edge" in str(error):
                continue
            else:
                raise error

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

    if not edge_path_index:
        return edge_fail_dictionary

    select_flows = flow_dataframe[flow_dataframe.index.isin(edge_path_index)]
    del edge_path_index

    first_edge_id = edge_failure_set[0]
    del edge_failure_set
    A = sorted(
        network_graph.clusters().subgraphs(),
        key=lambda l: len(l.es["edge_id"]),
        reverse=True,
    )
    access_flows = []
    for i in range(len(A)):
        network_graph = A[i]

        # TODO: most of the time this array probably doesn't change, we have one big island
        # but it's 9% of runtime, so perhaps create it once for the default condition,
        # and again if necessary
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
            access_flows, axis=0, sort=False, ignore_index=True
        )
        select_flows = pd.merge(
            select_flows,
            access_flows,
            how="left",
            on=["origin_id", "destination_id"],
        ).fillna(0)
    else:
        # TODO: SettingWithCopyWarning: 
        # A value is trying to be set on a copy of a slice from a DataFrame.
        # Try using .loc[row_indexer,col_indexer] = value instead
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


def read_flow_data(data_path: str):
    """
    Read combined flow data from disk ready for transport failure disruption.
    """

    logging.info("Reading trade sector list")
    with open(os.path.join(data_path, "trade_sectors.json"), "r") as fp:
        trade_sectors = json.load(fp)

    network_data: dict = {}
    flow_types = [f"trade_{sector}" for sector in trade_sectors] + ["labour"]
    for flow_type in flow_types:
        output_flow_dir = os.path.join(data_path, flow_type)
        logging.info(f"Reading {flow_type} network")
        network_df = gpd.read_parquet(os.path.join(output_flow_dir, "network.gpq"))
        network: ig.Graph = ig.Graph.TupleList(
            network_df.itertuples(index=False),
            edge_attrs=['edge_id', 'from_mode', 'to_mode', 'length_m', 'speed', 'time', 'geometry']
        )
        logging.info(f"Reading {flow_type} flows")
        flows = pd.read_parquet(os.path.join(output_flow_dir, "flows.pq"))
        logging.info(f"Reading {flow_type} edge indices")
        edge_indexes = pd.read_parquet(os.path.join(output_flow_dir, "edge_indexes.pq"))
        network_data[flow_type] = {
            "network": network,
            "flows": flows,
            "edge_indexes": edge_indexes.to_dict()["edge_indexes"],
        }

    logging.info("Reading combined flows")
    all_flows = pd.read_parquet(os.path.join(data_path, "all_flows.pq"))

    return network_data, all_flows, trade_sectors
