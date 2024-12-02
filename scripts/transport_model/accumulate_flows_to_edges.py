import logging
import multiprocessing
import sys

import geopandas as gpd
import pandas as pd


def init_worker(labour_flows_path: str):
    logging.info("Initialising worker")
    global od_pairs
    od_pairs = pd.read_parquet(labour_flows_path, columns=("edge_path", "working_trips", "GDP_to_trips"))
    return


def find_od_pairs_using_edge(edge_id) -> list[tuple[str, float, float]]:
    logging.info(edge_id)
    edge_trips = []
    for od_pair in od_pairs.itertuples():
        if edge_id in od_pair.edge_path:
            edge_trips.append(
                (
                    edge_id,
                    od_pair.working_trips,
                    od_pair.GDP_to_trips
                )
            )
    return edge_trips


def flatten(x: list[list]) -> list:
    return [item for sublist in x for item in sublist]


if __name__ == "__main__":

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

    _, n_cpu, network_path, labour_flows_path, edge_flows_path = sys.argv

    logging.info("Reading multi-modal network")
    network_edges = gpd.read_file(network_path, layer="edges")

    logging.info("Accumulating flows to edges")
    # ~50s per edge, 100k edges, therefore ~60d in serial!
    results = []
    edge_list = network_edges["edge_id"].tolist()
    with multiprocessing.Pool(
        processes=int(n_cpu),
        initializer=init_worker,
        initargs=(labour_flows_path,),
    ) as pool:
        results = pool.starmap(find_od_pairs_using_edge, [(edge_id,) for edge_id in edge_list])

    edge_trips = flatten(results)
    flows_by_edge = pd.DataFrame(
        edge_trips, columns=["edge_id", "working_trips", "GDP_to_trips"]
    )
    accumulated_flows = flows_by_edge.loc[:, ["edge_id", "working_trips", "GDP_to_trips"]].groupby("edge_id").sum().reset_index()

    logging.info("Writing accumulated flows to disk")
    accumulated_flows.to_parquet(edge_flows_path)