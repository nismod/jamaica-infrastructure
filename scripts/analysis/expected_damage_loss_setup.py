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
from analysis_utils import *
import subprocess 

#####################################
# READ MAIN DATA
#####################################

def main(config):
    processed_data_path = config['paths']['data']
    results_path = config['paths']['output']
    
    """Next we call the failure analysis script and loop through the falure scenarios
    """
    num_blocks = 14
    args = ["parallel",
            "-j", str(num_blocks),
            "--colsep", ",",
            "-a",
            "parameter_combinations.txt",
            "python",
            "expected_damages_losses_calculations_parallel.py",
            "{}"
            ]
    print (args)
    subprocess.run(args)
                                
if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)