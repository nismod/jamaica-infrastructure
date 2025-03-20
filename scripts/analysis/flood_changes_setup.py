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

    flood_thresholds = [1.0, 1.5, 2.0, 2.5]
    for ft in flood_thresholds:
        folder_name = f"flood_threshold_{str(ft).replace('.','p')}"
        results_folder = os.path.join(results_path, folder_name)
        if os.path.exists(results_folder) == False:
            os.mkdir(results_folder)

        damage_results_folder = f"{folder_name}/direct_damages"
        summary_folder = f"{folder_name}/direct_damages_summary"
        timeseries_results_folder = f"{folder_name}/loss_damage_timeseries"
        discounted_results_folder = f"{folder_name}/loss_damage_npvs"

        flood_protection_column = folder_name

        # Rework the networks file and write it to a new path
        networks = pd.read_csv(network_csv)
        networks = networks[networks["sector"] != "buildings"]
        network_csv = os.path.join(
            results_folder, "network_layers_hazard_intersections_details.csv"
        )
        networks.to_csv(network_csv, index=False)
        del networks

        # Rework the networks file and write it to a new path
        hazards = pd.read_csv(hazard_damage_parameters_csv)
        hazards = hazards[hazards["hazard_type"] == "flooding"]
        hazards["hazard_threshold"] = ft
        hazard_damage_parameters_csv = os.path.join(
            results_folder, "hazard_damage_parameters.csv.csv"
        )
        hazards.to_csv(hazard_damage_parameters_csv, index=False)
        del hazards

        parameter_combinations_file = "sensitivity_parameters.csv"

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
                        f"{damage_results_folder},{network_csv},{hazard_csv},{flood_protection_column},{pv[0]},{pv[1]},{pv[2]}\n"
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
