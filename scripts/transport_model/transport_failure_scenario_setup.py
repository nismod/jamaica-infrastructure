"""This script allows us to select and parallelise the failure analysis simulations on a server with multipe core processors
    The failure simulations selected here are for:
        Single node and edge failure of the transport networks 
    
    The script first reads in the single node/edge data which we want to fail 
        And then divides the numbers of selected scenarios into paritions which will be selected as initiating failure scenarios
        
        Example output - 0,90
        - First node to sample for initiating failure - The one at location 0 on the sample list
        - Last node to sample for initiating failure - The one at location 90 on the sample list  
        
        We are mainly selecting the first 90 nodes of the road and rail network one-by-one and failing them
        
        All partitions are first written into text files named as parallel_network_scenario_selection.txt 
        Example: parallel_transport_scenario_selection.txt 

        Example output in file: parallel_transport_scenario_selection.txt 
            0,90
            90,181
            181,272
            272,363
            363,453
            453,544
            544,635
        
        Each of these lines is a batch of scenarios that are run on different processors in parallel
"""

import datetime
import os
import pandas as pd
import subprocess
import sys

import geopandas as gpd
import numpy as np
from tqdm import tqdm

from jamaica_infrastructure.transport.utils import load_config

#####################################
# READ MAIN DATA
#####################################


def main(config, n_cpu):
    processed_data_path = config["paths"]["data"]
    results_path = config["paths"]["output"]

    num_blocks = 200  # Number of partitions of the networks nodes created for parallel processing
    edges = gpd.read_file(
        os.path.join(
            processed_data_path, "networks", "transport", "multi_modal_network.gpkg"
        ),
        layer="edges",
    )
    rail_edges = edges[(edges["from_mode"] == "rail") & (edges["to_mode"] == "rail")][
        "edge_id"
    ].values.tolist()
    road_edges = edges[(edges["from_mode"] == "road") & (edges["to_mode"] == "road")][
        "edge_id"
    ].values.tolist()
    edge_fail = rail_edges + road_edges

    num_values = np.linspace(0, len(edge_fail) - 1, num_blocks)
    scenarios_path = os.path.join(results_path, "transport_failures", "parallel_transport_scenario_selection.txt")
    with open(scenarios_path, "w+") as f:
        for n in range(len(num_values) - 1):
            f.write(
                "{},{},{}\n".format(
                    "fail edges", int(num_values[n]), int(num_values[n + 1])
                )
            )

    scenarios_resampled_path = os.path.join(
        results_path,
        "transport_failures",
        "parallel_transport_scenario_selection_resample.txt"
    )
    with open(scenarios_resampled_path, "w+") as f:
        with open(scenarios_path) as t:
            for line in t:
                line = line.strip("\n")
                file = os.path.join(
                    results_path,
                    "transport_failures",
                    f"single_link_failures_scenarios_{line.split(',')[1]}_{line.split(',')[2]}.csv",
                )
                print(file)
                if os.path.exists(file) is False:
                    f.write(f"{line}\n")

    """Next we call the failure analysis script and loop through the failure scenarios
    """
    args = [
        "parallel",
        "-j",
        str(n_cpu),
        "--colsep",
        ",",
        "-a",
        scenarios_path,
        "python",
        "transport_failure_analysis.py",
        "{}",
    ]
    print(args)
    subprocess.run(args)

    # signal to snakemake that the job is complete
    with open(os.path.join(results_path, "transport_failures", "transport_failures.flag"), "w") as fp:
        fp.write(datetime.datetime.now())


if __name__ == "__main__":

    _, n_cpu = sys.argv

    CONFIG = load_config()
    main(CONFIG, int(n_cpu))