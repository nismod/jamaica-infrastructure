import logging
import multiprocessing
import sys

import geopandas as gpd
import pandas as pd
from tqdm import tqdm


def increment_edge_flows(edge_flows: pd.DataFrame, od_pairs_chunk: pd.DataFrame) -> pd.DataFrame:
    for od_pair in tqdm(od_pairs_chunk.itertuples(), total=len(od_pairs_chunk)):
        edge_flows.loc[od_pair.edge_path, "working_trips"] += od_pair.working_trips
        edge_flows.loc[od_pair.edge_path, "GDP_to_trips"] += od_pair.GDP_to_trips
    return edge_flows


if __name__ == "__main__":

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

    _, n_cpu, network_path, labour_flows_path, edge_flows_path, edge_flows_geometry_path = sys.argv
    n_cpu = int(n_cpu)

    logging.info("Reading multi-modal network")
    network_edges = gpd.read_file(network_path, layer="edges")

    logging.info("Reading flows")
    od_pairs = pd.read_parquet(labour_flows_path, columns=("edge_path", "working_trips", "GDP_to_trips"))

    # empty edge flows dataframe to accumulate flows to
    edge_flows = network_edges.loc[:, ["edge_id"]].copy()
    edge_flows["working_trips"] = 0.0
    edge_flows["GDP_to_trips"] = 0.0
    edge_flows = edge_flows.set_index("edge_id")

    chunk_size = int(len(od_pairs) // n_cpu + 1)
    logging.info(f"Splitting OD pairs into {n_cpu} chunks of ~{chunk_size}")
    args = []
    for i in range(n_cpu):
        args.append((edge_flows.copy(), od_pairs.iloc[i * chunk_size: (i + 1) * chunk_size].copy()))

    logging.info("Accumulating flows to edges")
    with multiprocessing.Pool(processes=n_cpu) as pool:
        results = pool.starmap(increment_edge_flows, args)
    
    # combine the chunks, summing across edges
    edge_flows, *other_edge_flows = results
    for other in other_edge_flows:
        edge_flows += other

    logging.info("Writing accumulated flows to disk")
    edge_flows = edge_flows.reset_index()
    edge_flows.to_parquet(edge_flows_path)
    # join with geometry
    edge_flows_geometry = gpd.GeoDataFrame(
        edge_flows.merge(network_edges.loc[:, ["edge_id", "geometry"]], on="edge_id")
    )
    edge_flows_geometry.to_parquet(edge_flows_geometry_path)
