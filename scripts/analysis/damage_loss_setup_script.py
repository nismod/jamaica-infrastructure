"""This script allows us to select and parallelise the Damage and Loss estimations on a server with multiple core processors
"""

import os
from SALib.sample import morris
from analysis_utils import *
import subprocess


def main(config, network_csv, hazard_csv, n_cpu, sensitivity_analysis):
    processed_data_path = config["paths"]["data"]
    results_path = config["paths"]["output"]

    network_csv = os.path.join(
        processed_data_path,
        "networks",
        "network_layers_hazard_intersections_details.csv",
    )
    damage_curves_csv = os.path.join(
        processed_data_path, "damage_curves", "asset_damage_curve_mapping.csv"
    )
    hazard_damage_parameters_csv = os.path.join(
        processed_data_path, "damage_curves", "hazard_damage_parameters.csv"
    )

    damage_results_folder = os.path.join(results_path, "direct_damages")
    summary_folder = os.path.join(results_path, "direct_damages_summary")
    timeseries_results_folder = os.path.join(results_path, "loss_damage_timeseries")
    discounted_results_folder = os.path.join(results_path, "loss_damage_npvs")
    for folder in (damage_results_folder, summary_folder, timeseries_results_folder, discounted_results_folder):
        if not os.path.exists(folder):
            os.makedirs(folder)

    flood_protection_column = "None"

    parameter_combinations_file = os.path.join(
        processed_data_path, "parameter_combinations.txt"
    )

    if sensitivity_analysis:
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

    else:
        param_values = [(1.0, 1.0)]  # one sample, no variation

    with open(parameter_combinations_file, "w+") as f:
        for p in range(len(param_values)):
            f.write(f"{p},{param_values[p][0]},{param_values[p][1]}\n")

    damage_results_fname = os.path.join(damage_results_folder, "damage_results.txt")
    with open(damage_results_fname, "w+") as f:
        with open(parameter_combinations_file, "r") as r:
            for p in r:
                pv = p.split(",")
                f.write(
                    f"{damage_results_folder},{network_csv},{hazard_csv},{damage_curves_csv},{hazard_damage_parameters_csv},{pv[0]},{pv[1]},{pv[2]}\n"
                )

    """Call the failure analysis script and loop through the failure scenarios"""
    args = [
        "parallel",
        "--halt",
        "now,fail=1",
        "--lb",
        "-j",
        str(n_cpu),
        "--colsep",
        ",",
        "-a",
        damage_results_fname,
        "python",
        "scripts/analysis/damage_calculations.py",
        "{}",
    ]
    print("* Start the processing of damage calculations")
    print(args)
    subprocess.check_output(args)

    ead_eael_results_fname = os.path.join(damage_results_folder, "ead_eael_results.txt")
    with open(ead_eael_results_fname, "w+") as f:
        with open(parameter_combinations_file, "r") as r:
            for p in r:
                pv = p.split(",")
                f.write(
                    f"{damage_results_folder},{network_csv},{hazard_csv},{flood_protection_column},{pv[0]},{pv[1]},{pv[2]}"
                )


if __name__ == "__main__":
    CONFIG = load_config()
    network_csv, hazard_csv, n_cpu, sensitivity_analysis = sys.argv[1:]
    main(CONFIG, network_csv, hazard_csv, int(n_cpu), bool(int(sensitivity_analysis)))
