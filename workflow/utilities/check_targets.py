"""
This script checks to see if there's a sequence of rules available for creating
the chain of rule inputs and outputs which ultimately may generate a desired
target file.

It does not try and generate any of the intermediate or target files.

Run from the root of the repository:
$ python workflow/utilities/check_targets.py

If .tmp already exists, it will need to be deleted first.
"""

import os
import subprocess
import re
import sys


EXPOSURE = [
    "results/direct_damages_summary_uids/airport_polygon_areas_exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_2.6__epoch_2050__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_2.6__epoch_2100__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_4.5__epoch_2030__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_4.5__epoch_2050__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_4.5__epoch_2070__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_4.5__epoch_2100__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_8.5__epoch_2030__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_8.5__epoch_2050__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_8.5__epoch_2070__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_8.5__epoch_2100__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_baseline__epoch_2010__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_cyclone__rcp_4.5__epoch_2050__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_cyclone__rcp_4.5__epoch_2100__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_cyclone__rcp_8.5__epoch_2050__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_cyclone__rcp_8.5__epoch_2100__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_cyclone__rcp_baseline__epoch_2010__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_fluvial__rcp_2.6__epoch_2050__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_fluvial__rcp_2.6__epoch_2080__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_fluvial__rcp_4.5__epoch_2050__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_fluvial__rcp_4.5__epoch_2080__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_fluvial__rcp_8.5__epoch_2050__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_fluvial__rcp_8.5__epoch_2080__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_fluvial__rcp_baseline__epoch_2010__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_surface__rcp_2.6__epoch_2050__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_surface__rcp_2.6__epoch_2080__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_surface__rcp_4.5__epoch_2050__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_surface__rcp_4.5__epoch_2080__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_surface__rcp_8.5__epoch_2050__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_surface__rcp_8.5__epoch_2080__exposures.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_surface__rcp_baseline__epoch_2010__exposures.parquet",
    "results/direct_damages_summary_uids/electricity_network_v3.1_edges_exposures.parquet",
    "results/direct_damages_summary_uids/electricity_network_v3.1_nodes_exposures.parquet",
    "results/direct_damages_summary_uids/irrigation_assets_NIC_edges_exposures.parquet",
    "results/direct_damages_summary_uids/irrigation_assets_NIC_nodes_exposures.parquet",
    "results/direct_damages_summary_uids/pipelines_NWC_edges_exposures.parquet",
    "results/direct_damages_summary_uids/port_polygon_areas_exposures.parquet",
    "results/direct_damages_summary_uids/potable_facilities_NWC_nodes_exposures.parquet",
    "results/direct_damages_summary_uids/rail_edges_exposures.parquet",
    "results/direct_damages_summary_uids/rail_nodes_exposures.parquet",
    "results/direct_damages_summary_uids/roads_edges_exposures.parquet",
    "results/direct_damages_summary_uids/roads_nodes_exposures.parquet",
    "results/direct_damages_summary_uids/waste_water_facilities_NWC_nodes_exposures.parquet",
]

DAMAGES = [
    "results/direct_damages_summary_uids/airport_polygon_areas_damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_2.6__epoch_2050__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_2.6__epoch_2100__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_4.5__epoch_2030__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_4.5__epoch_2050__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_4.5__epoch_2070__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_4.5__epoch_2100__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_8.5__epoch_2030__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_8.5__epoch_2050__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_8.5__epoch_2070__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_8.5__epoch_2100__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_baseline__epoch_2010__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_cyclone__rcp_4.5__epoch_2050__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_cyclone__rcp_4.5__epoch_2100__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_cyclone__rcp_8.5__epoch_2050__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_cyclone__rcp_8.5__epoch_2100__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_cyclone__rcp_baseline__epoch_2010__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_fluvial__rcp_2.6__epoch_2050__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_fluvial__rcp_2.6__epoch_2080__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_fluvial__rcp_4.5__epoch_2050__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_fluvial__rcp_4.5__epoch_2080__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_fluvial__rcp_8.5__epoch_2050__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_fluvial__rcp_8.5__epoch_2080__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_fluvial__rcp_baseline__epoch_2010__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_surface__rcp_2.6__epoch_2050__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_surface__rcp_2.6__epoch_2080__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_surface__rcp_4.5__epoch_2050__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_surface__rcp_4.5__epoch_2080__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_surface__rcp_8.5__epoch_2050__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_surface__rcp_8.5__epoch_2080__damages.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_surface__rcp_baseline__epoch_2010__damages.parquet",
    "results/direct_damages_summary_uids/electricity_network_v3.1_edges_damages.parquet",
    "results/direct_damages_summary_uids/electricity_network_v3.1_nodes_damages.parquet",
    "results/direct_damages_summary_uids/irrigation_assets_NIC_edges_damages.parquet",
    "results/direct_damages_summary_uids/irrigation_assets_NIC_nodes_damages.parquet",
    "results/direct_damages_summary_uids/pipelines_NWC_edges_damages.parquet",
    "results/direct_damages_summary_uids/port_polygon_areas_damages.parquet",
    "results/direct_damages_summary_uids/potable_facilities_NWC_nodes_damages.parquet",
    "results/direct_damages_summary_uids/rail_edges_damages.parquet",
    "results/direct_damages_summary_uids/rail_nodes_damages.parquet",
    "results/direct_damages_summary_uids/roads_edges_damages.parquet",
    "results/direct_damages_summary_uids/roads_nodes_damages.parquet",
    "results/direct_damages_summary_uids/waste_water_facilities_NWC_nodes_damages.parquet",
]

LOSSES = [
    "results/direct_damages_summary_uids/airport_polygon_areas_losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_2.6__epoch_2050__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_2.6__epoch_2100__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_4.5__epoch_2030__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_4.5__epoch_2050__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_4.5__epoch_2070__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_4.5__epoch_2100__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_8.5__epoch_2030__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_8.5__epoch_2050__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_8.5__epoch_2070__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_8.5__epoch_2100__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_coastal__rcp_baseline__epoch_2010__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_cyclone__rcp_4.5__epoch_2050__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_cyclone__rcp_4.5__epoch_2100__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_cyclone__rcp_8.5__epoch_2050__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_cyclone__rcp_8.5__epoch_2100__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_cyclone__rcp_baseline__epoch_2010__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_fluvial__rcp_2.6__epoch_2050__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_fluvial__rcp_2.6__epoch_2080__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_fluvial__rcp_4.5__epoch_2050__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_fluvial__rcp_4.5__epoch_2080__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_fluvial__rcp_8.5__epoch_2050__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_fluvial__rcp_8.5__epoch_2080__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_fluvial__rcp_baseline__epoch_2010__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_surface__rcp_2.6__epoch_2050__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_surface__rcp_2.6__epoch_2080__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_surface__rcp_4.5__epoch_2050__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_surface__rcp_4.5__epoch_2080__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_surface__rcp_8.5__epoch_2050__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_surface__rcp_8.5__epoch_2080__losses.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_surface__rcp_baseline__epoch_2010__losses.parquet",
    "results/direct_damages_summary_uids/electricity_network_v3.1_edges_losses.parquet",
    "results/direct_damages_summary_uids/electricity_network_v3.1_nodes_losses.parquet",
    "results/direct_damages_summary_uids/irrigation_assets_NIC_edges_losses.parquet",
    "results/direct_damages_summary_uids/irrigation_assets_NIC_nodes_losses.parquet",
    "results/direct_damages_summary_uids/pipelines_NWC_edges_losses.parquet",
    "results/direct_damages_summary_uids/port_polygon_areas_losses.parquet",
    "results/direct_damages_summary_uids/potable_facilities_NWC_nodes_losses.parquet",
    "results/direct_damages_summary_uids/rail_edges_losses.parquet",
    "results/direct_damages_summary_uids/rail_nodes_losses.parquet",
    "results/direct_damages_summary_uids/roads_edges_losses.parquet",
    "results/direct_damages_summary_uids/roads_nodes_losses.parquet",
    "results/direct_damages_summary_uids/waste_water_facilities_NWC_nodes_losses.parquet",
]

EAD_EAEL = [
    "results/direct_damages_summary_uids/airport_polygon_areas_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/electricity_network_v3.1_edges_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/electricity_network_v3.1_nodes_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/irrigation_assets_NIC_edges_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/irrigation_assets_NIC_nodes_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/pipelines_NWC_edges_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/port_polygon_areas_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/potable_facilities_NWC_nodes_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/rail_edges_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/rail_nodes_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/roads_edges_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/roads_nodes_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/waste_water_facilities_NWC_nodes_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/airport_polygon_areas_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/electricity_network_v3.1_edges_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/electricity_network_v3.1_nodes_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/irrigation_assets_NIC_edges_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/irrigation_assets_NIC_nodes_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/pipelines_NWC_edges_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/port_polygon_areas_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/potable_facilities_NWC_nodes_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/rail_edges_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/rail_nodes_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/roads_edges_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/roads_nodes_EAD_EAEL.parquet",
    "results/direct_damages_summary_uids/waste_water_facilities_NWC_nodes_EAD_EAEL.parquet",
]

ADAPTATIONS = [
    'results/adaptation_benefits_costs_bcr/flooding_airport_polygon_areas_adaptation_benefits_costs_bcr.csv',
    'results/adaptation_benefits_costs_bcr/flooding_airport_polygon_areas_adaptation_costs_avoided_EAD_EAEL.csv',
    'results/adaptation_benefits_costs_bcr/flooding_electricity_network_v3.1_nodes_adaptation_benefits_costs_bcr.csv',
    'results/adaptation_benefits_costs_bcr/flooding_electricity_network_v3.1_edges_adaptation_costs_avoided_EAD_EAEL.csv',
    'results/adaptation_benefits_costs_bcr/flooding_irrigation_assets_NIC_nodes_adaptation_benefits_costs_bcr.csv',
    'results/adaptation_benefits_costs_bcr/flooding_irrigation_assets_NIC_nodes_adaptation_costs_avoided_EAD_EAEL.csv',
    'results/adaptation_benefits_costs_bcr/flooding_port_polygon_areas_adaptation_benefits_costs_bcr.csv',
    'results/adaptation_benefits_costs_bcr/flooding_port_polygon_areas_adaptation_costs_avoided_EAD_EAEL.csv',
    'results/adaptation_benefits_costs_bcr/flooding_potable_facilities_NWC_nodes_adaptation_benefits_costs_bcr.csv',
    'results/adaptation_benefits_costs_bcr/flooding_potable_facilities_NWC_nodes_adaptation_costs_avoided_EAD_EAEL.csv',
    'results/adaptation_benefits_costs_bcr/flooding_rail_edges_adaptation_benefits_costs_bcr.csv',
    'results/adaptation_benefits_costs_bcr/flooding_rail_edges_adaptation_costs_avoided_EAD_EAEL.csv',
    'results/adaptation_benefits_costs_bcr/flooding_roads_edges_adaptation_benefits_costs_bcr.csv',
    'results/adaptation_benefits_costs_bcr/flooding_roads_edges_adaptation_costs_avoided_EAD_EAEL.csv',
    'results/adaptation_benefits_costs_bcr/flooding_roads_nodes_adaptation_benefits_costs_bcr.csv',
    'results/adaptation_benefits_costs_bcr/flooding_roads_nodes_adaptation_costs_avoided_EAD_EAEL.csv',
    'results/adaptation_benefits_costs_bcr/flooding_waste_water_facilities_NWC_nodes_adaptation_benefits_costs_bcr.csv',
    'results/adaptation_benefits_costs_bcr/flooding_waste_water_facilities_NWC_nodes_adaptation_costs_avoided_EAD_EAEL.csv',
    'results/adaptation_benefits_costs_bcr/TC_electricity_network_v3.1_nodes_adaptation_benefits_costs_bcr.csv',
    'results/adaptation_benefits_costs_bcr/TC_electricity_network_v3.1_edges_adaptation_costs_avoided_EAD_EAEL.csv',
    'results/adaptation_benefits_costs_bcr/coastal_rail_nodes_adaptation_benefits_costs_bcr.csv',
    'results/adaptation_benefits_costs_bcr/coastal_rail_nodes_adaptation_costs_avoided_EAD_EAEL.csv',
    'results/adaptation_benefits_costs_bcr/coastal_rail_edges_adaptation_benefits_costs_bcr.csv',
    'results/adaptation_benefits_costs_bcr/coastal_rail_edges_adaptation_costs_avoided_EAD_EAEL.csv',
    'results/adaptation_benefits_costs_bcr/coastal_airport_polygon_areas_adaptation_benefits_costs_bcr.csv',
    'results/adaptation_benefits_costs_bcr/coastal_airport_polygon_areas_adaptation_costs_avoided_EAD_EAEL.csv',
    'results/adaptation_benefits_costs_bcr/coastal_port_polygon_areas_adaptation_benefits_costs_bcr.csv',
    'results/adaptation_benefits_costs_bcr/coastal_port_polygon_areas_adaptation_costs_avoided_EAD_EAEL.csv',
    'results/adaptation_benefits_costs_bcr/coastal_irrigation_assets_NIC_nodes_adaptation_benefits_costs_bcr.csv',
    'results/adaptation_benefits_costs_bcr/coastal_irrigation_assets_NIC_nodes_adaptation_costs_avoided_EAD_EAEL.csv',
    'results/adaptation_benefits_costs_bcr/coastal_irrigation_assets_NIC_edges_adaptation_benefits_costs_bcr.csv',
    'results/adaptation_benefits_costs_bcr/coastal_irrigation_assets_NIC_edges_adaptation_costs_avoided_EAD_EAEL.csv',
    'results/adaptation_benefits_costs_bcr/coastal_potable_facilities_NWC_nodes_adaptation_benefits_costs_bcr.csv',
    'results/adaptation_benefits_costs_bcr/coastal_potable_facilities_NWC_nodes_adaptation_costs_avoided_EAD_EAEL.csv',
    'results/adaptation_benefits_costs_bcr/coastal_waste_water_facilities_NWC_nodes_adaptation_benefits_costs_bcr.csv',
    'results/adaptation_benefits_costs_bcr/coastal_waste_water_facilities_NWC_nodes_adaptation_costs_avoided_EAD_EAEL.csv',
    'results/adaptation_benefits_costs_bcr/coastal_pipelines_NWC_edges_adaptation_benefits_costs_bcr.csv',
    'results/adaptation_benefits_costs_bcr/coastal_pipelines_NWC_edges_adaptation_costs_avoided_EAD_EAEL.csv',
    'results/adaptation_benefits_costs_bcr/coastal_roads_nodes_adaptation_benefits_costs_bcr.csv',
    'results/adaptation_benefits_costs_bcr/coastal_roads_nodes_adaptation_costs_avoided_EAD_EAEL.csv',
    'results/adaptation_benefits_costs_bcr/coastal_roads_edges_adaptation_benefits_costs_bcr.csv',
    'results/adaptation_benefits_costs_bcr/coastal_roads_edges_adaptation_costs_avoided_EAD_EAEL.csv',
    'results/adaptation_benefits_costs_bcr/coastal_electricity_network_v3.1_nodes_adaptation_benefits_costs_bcr.csv',
    'results/adaptation_benefits_costs_bcr/coastal_electricity_network_v3.1_nodes_adaptation_costs_avoided_EAD_EAEL.csv',
    'results/adaptation_benefits_costs_bcr/coastal_electricity_network_v3.1_edges_adaptation_benefits_costs_bcr.csv',
    'results/adaptation_benefits_costs_bcr/coastal_electricity_network_v3.1_edges_adaptation_costs_avoided_EAD_EAEL.csv',
    'results/adaptation_benefits_costs_bcr/coastal_buildings_assigned_economic_activity_areas_adaptation_benefits_costs_bcr.csv',
    'results/adaptation_benefits_costs_bcr/coastal_buildings_assigned_economic_activity_areas_adaptation_costs_avoided_EAD_EAEL.csv',
]

HOTSPOTS = [
    "results/hotspots/exposure/all_sectors.tiff",
    "results/hotspots/exposure/water.tiff",
    "results/hotspots/exposure/transport.tiff",
    "results/hotspots/exposure/energy.tiff",
    "results/hotspots/transport/isolation_loss_smoothed.tiff",
    "results/hotspots/transport/rerouting_loss_smoothed.tiff",
]

TARGETS = EXPOSURE + DAMAGES + LOSSES + EAD_EAEL + ADAPTATIONS + HOTSPOTS

if __name__ == "__main__":
    n_okay = 0
    n_warn = 0

    if os.path.exists(".tmp"):
        sys.exit("Please remove .tmp folder before running.")
    else:
        os.makedirs(".tmp/rulegraphs")

    for target in TARGETS:

        # run snakemake to check the rulegraph
        # print the output if it fails
        print(f"Checking {target}... ", end="")
        run = subprocess.run(
            f"snakemake --cores 1 --rulegraph {target}",
            capture_output=True,
            shell=True,
        )

        if run.returncode == 0:
            n_okay += 1
            matches = re.findall("(Warning: .+)", run.stdout.decode(), re.IGNORECASE)
            if len(matches) > 0:
                print(f"{len(matches)} warnings:")
                print("\n".join(matches))
                n_warn += 1
            else:
                fname = os.path.basename(target).replace('.parquet', '.png').replace('.csv', '.png')

                # plot the rulegraph and write to disk
                subprocess.run(
                    ["dot", "-Tpng", "-o", f".tmp/rulegraphs/{fname}"],
                    check=True,
                    input=run.stdout,
                )

                print("OK")
        else:
            print("FAILED:")
            print(run.stdout.decode())
            print(run.stderr.decode())

    if n_okay == len(TARGETS):
        print("All targets pass dry-run test.")
    else:
        print(f"{n_okay}/{len(TARGETS)} targets pass dry-run test.")
