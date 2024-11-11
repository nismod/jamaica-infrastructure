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
import sys
import ujson
from SALib.sample import morris
import SALib.analyze.morris
from analysis_utils import *
import subprocess

#####################################
# READ MAIN DATA
#####################################


def main(config):
    processed_data_path = config["paths"]["data"]
    results_path = config["paths"]["output"]

    # Set up problem for sensitivity analysis
    problem = {
        "num_vars": 2,
        "names": ["cost_uncertainty_parameter", "damage_uncertainty_parameter"],
        "bounds": [[0.0, 1.0], [0.0, 1.0]],
    }

    # And create parameter values
    param_values = morris.sample(
        problem, 10, num_levels=4, optimal_trajectories=8, local_optimization=False
    )
    param_values = list(set([(p[0], p[1]) for p in param_values]))
    with open("parameter_combinations.txt", "w+") as f:
        # f.write("parameter_set cost_uncertainty_parameter damage_uncertainty_parameter\n")
        for p in range(len(param_values)):
            f.write(f"{p},{param_values[p][0]},{param_values[p][1]}\n")

    f.close()
    num_blocks = len(param_values)
    """Next we call the failure analysis script and loop through the falure scenarios
    """
    args = [
        "parallel",
        "-j",
        str(num_blocks),
        "--colsep",
        ",",
        "-a",
        "parameter_combinations.txt",
        "python",
        "direct_damage_calculations_parallel.py",
        "{}",
    ]
    print(args)
    subprocess.run(args)


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
