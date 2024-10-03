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

import os
import subprocess

import pandas as pd
import geopandas as gpd
from tqdm import tqdm

from jamaica_infrastructure.transport.utils import load_config, np

#####################################
# READ MAIN DATA
#####################################


def main(config):
    processed_data_path = config["paths"]["data"]
    results_path = config["paths"]["output"]

    num_blocks = (
        20  # Number of partitions of the networks nodes created for parallel processing
    )
    bridges = gpd.read_file(
        os.path.join(processed_data_path, "networks", "transport", "roads.gpkg"),
        layer="nodes",
    )
    bridges = bridges[bridges["asset_type"] == "bridge"]["node_id"].values.tolist()

    num_values = np.linspace(0, len(bridges) - 1, num_blocks)
    with open("parallel_bridge_scenario_selection.txt", "w+") as f:
        for n in range(len(num_values) - 1):
            f.write(
                "{},{},{}\n".format(
                    "fail bridges", int(num_values[n]), int(num_values[n + 1])
                )
            )

    f.close()

    # with open("parallel_transport_scenario_selection_resample.txt","w+") as f:
    #     with open("parallel_transport_scenario_selection.txt") as t:
    #         for line in t:
    #             line = line.strip('\n')
    #             file = os.path.join(results_path,
    #                             "transport_failures",
    #                             f"single_link_failures_scenarios_{line.split(',')[1]}_{line.split(',')[2]}.csv")
    #             print (file)
    #             if os.path.exists(file) is False:
    #                 f.write(f"{line}\n")

    # f.close()

    """Next we call the failure analysis script and loop through the falure scenarios
    """
    args = [
        "parallel",
        "-j",
        str(num_blocks),
        "--colsep",
        ",",
        "-a",
        "parallel_bridge_scenario_selection.txt",
        "python",
        "roads_bridges_failure_analysis.py",
        "{}",
    ]
    print(args)
    subprocess.run(args)


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
