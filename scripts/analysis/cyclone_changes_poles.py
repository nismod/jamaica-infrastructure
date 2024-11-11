"""This script allows us to select and parallelise the Damage and Loss estimations on a server with multiple core processors
"""

import os
import sys
import ujson
from SALib.sample import morris
import SALib.analyze.morris
from analysis_utils import *
import subprocess


def main(config):
    processed_data_path = config["paths"]["data"]
    results_path = config["paths"]["output"]

    network_csv = os.path.join(
        processed_data_path,
        "networks",
        "network_layers_hazard_intersections_details.csv",
    )
    hazard_csv = os.path.join(processed_data_path, "hazards", "hazard_layers.csv")
    damage_curves_csv = os.path.join(
        processed_data_path, "damage_curves", "asset_damage_curve_mapping.csv"
    )
    hazard_damage_parameters_csv = os.path.join(
        processed_data_path, "damage_curves", "hazard_damage_parameters.csv"
    )

    cyclone_damage_curve_change = [0.76]
    for ft in cyclone_damage_curve_change:
        folder_name = f"cyclone_damage_curve_change_{str(ft).replace('.','p')}"
        results_folder = os.path.join(results_path, folder_name)
        if os.path.exists(results_folder) == False:
            os.mkdir(results_folder)

        damage_results_folder = f"{folder_name}/direct_damages"
        summary_folder = f"{folder_name}/direct_damages_summary"
        timeseries_results_folder = f"{folder_name}/loss_damage_timeseries"
        discounted_results_folder = f"{folder_name}/loss_damage_npvs"

        cyclone_protection_column = folder_name

        # Rework the networks file and write it to a new path
        networks = pd.read_csv(network_csv)
        networks = networks[
            (networks["sector"] == "electricity") & (networks["asset_layer"] == "nodes")
        ]
        network_csv = os.path.join(
            results_folder, "network_layers_hazard_intersections_details.csv"
        )
        networks.to_csv(network_csv, index=False)
        del networks

        # Rework the networks file and write it to a new path
        hazards = pd.read_csv(hazard_damage_parameters_csv)
        hazards = hazards[hazards["hazard_type"] == "cyclone"]
        hazards["uplift_factor"] = -1.0 * ft
        hazard_damage_parameters_csv = os.path.join(
            results_folder, "hazard_damage_parameters.csv"
        )
        hazards.to_csv(hazard_damage_parameters_csv, index=False)
        del hazards

        # Rework the damage curves file
        damage_curves = pd.read_csv(damage_curves_csv)
        damage_curves = damage_curves[damage_curves["asset_sheet"] == "poles"]
        damage_curves_csv = os.path.join(
            results_folder, "asset_damage_curve_mapping.csv"
        )
        damage_curves.to_csv(damage_curves_csv, index=False)
        del damage_curves

        parameter_combinations_file = "parameter_combinations.txt"
        generate_new_parameters = False
        if generate_new_parameters is True:
            # Set up problem for sensitivity analysis
            problem = {
                "num_vars": 2,
                "names": ["cost_uncertainty_parameter", "damage_uncertainty_parameter"],
                "bounds": [[0.0, 1.0], [0.0, 1.0]],
            }

            # And create parameter values
            param_values = morris.sample(
                problem,
                10,
                num_levels=4,
                optimal_trajectories=8,
                local_optimization=False,
            )
            param_values = list(set([(p[0], p[1]) for p in param_values]))
            with open(parameter_combinations_file, "w+") as f:
                for p in range(len(param_values)):
                    f.write(f"{p},{param_values[p][0]},{param_values[p][1]}\n")

            f.close()

        with open("damage_results.txt", "w+") as f:
            with open(parameter_combinations_file, "r") as r:
                for p in r:
                    pv = p.split(",")
                    f.write(
                        f"{damage_results_folder},{network_csv},{hazard_csv},{damage_curves_csv},{hazard_damage_parameters_csv},{pv[0]},{pv[1]},{pv[2]}\n"
                    )

        f.close()

        num_blocks = len(param_values)
        """Next we call the failure analysis script and loop through the failure scenarios
        """
        args = [
            "parallel",
            "-j",
            str(num_blocks),
            "--colsep",
            ",",
            "-a",
            "damage_results.txt",
            "python",
            "damage_calculations.py",
            "{}",
        ]
        print("* Start the processing of damage calculations")
        print(args)
        subprocess.run(args)

        with open("ead_eael_results.txt", "w+") as f:
            with open(parameter_combinations_file, "r") as r:
                for p in r:
                    pv = p.split(",")
                    f.write(
                        f"{damage_results_folder},{network_csv},{hazard_csv},{cyclone_protection_column},{pv[0]},{pv[1]},{pv[2]}\n"
                    )

        f.close()

        """Next we call the EAD and EAEL analysis script and loop through the failure results
        """
        args = [
            "parallel",
            "-j",
            str(num_blocks),
            "--colsep",
            ",",
            "-a",
            "ead_eael_results.txt",
            "python",
            "ead_eael_calculations.py",
            "{}",
        ]
        print("* Start the processing of EAD and EAEL calculations")
        print(args)
        subprocess.run(args)

        """Next we call the summary scripts
        """
        args = [
            "python",
            "damage_loss_summarised.py",
            f"{damage_results_folder}",
            f"{summary_folder}",
            f"{network_csv}",
            f"{parameter_combinations_file}",
        ]
        print("* Start the processing of summarising damage results")
        print(args)
        subprocess.run(args)

        """Next we call the timeseries and NPV scripts
        """
        args = [
            "python",
            "damage_loss_timeseries_and_npv.py",
            f"{summary_folder}",
            f"{timeseries_results_folder}",
            f"{discounted_results_folder}",
            f"{network_csv}",
        ]
        print("* Start the processing of timeseries and NPV calculations")
        print(args)
        subprocess.run(args)


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
