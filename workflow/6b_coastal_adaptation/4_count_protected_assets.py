import logging
import os

import click
import geopandas as gpd
import pandas as pd


def process_network_assets(network_info, flood_polygons, data_path):
    """Process a single network and count assets and their costs intersecting with flood polygons."""
    fname = os.path.join(data_path, network_info["path"])
    layer_type = network_info["asset_layer"]
    asset = network_info["asset_gpkg"]

    def is_empty(val):
        return val is None or val == "" or str(val).lower() == "none" or pd.isna(val)

    mean_cost_col = network_info.get("asset_mean_cost_column")
    if is_empty(mean_cost_col):
        mean_cost_col = network_info.get("asset_max_cost_column")
    if is_empty(mean_cost_col):
        mean_cost_col = network_info.get("asset_min_cost_column")

    ref = f"{asset}_{layer_type}"

    # Initialize results dictionaries
    count_results = {str(pid): 0 for pid in flood_polygons["id"]}
    cost_results = {str(pid): 0.0 for pid in flood_polygons["id"]}

    logging.info(f"Processing network '{ref}'")

    try:
        # Load the asset layer
        assets = gpd.read_file(fname, layer=layer_type)
        if assets.empty:
            logging.warning(f"No assets found in layer {layer_type}")
            return count_results, cost_results

        logging.info(f"Found {len(assets)} assets in network '{ref}'")

        # Fix invalid geometries
        assets["geometry"] = assets["geometry"].apply(
            lambda geom: geom if geom.is_valid else geom.buffer(0)
        )
        flood_copy = flood_polygons.copy()
        flood_copy["geometry"] = flood_copy["geometry"].apply(
            lambda geom: geom if geom.is_valid else geom.buffer(0)
        )

        # Process each polygon
        for _, poly in flood_copy.iterrows():
            polygon_id = str(poly["id"])
            intersecting_assets = assets[assets.geometry.intersects(poly["geometry"])]

            count = len(intersecting_assets)
            count_results[polygon_id] = count

            # Sum costs from the mean cost column
            if mean_cost_col in intersecting_assets.columns:
                cost = intersecting_assets[mean_cost_col].fillna(0).sum()
            else:
                logging.warning(f"Column '{mean_cost_col}' not found in layer '{ref}'")
                cost = 0.0

            cost_results[polygon_id] = cost

            if count > 0:
                logging.debug(
                    f"Polygon {polygon_id}: {count} assets, total cost = {cost}"
                )

    except Exception as e:
        logging.error(f"Error processing network {ref}: {str(e)}")
        # logging.error(traceback.format_exc())

    total_count = sum(count_results.values())
    total_cost = sum(cost_results.values())
    logging.info(
        f"Network '{ref}': {total_count} total assets, {total_cost} total cost"
    )

    return count_results, cost_results


def count_assets(RCP, RP, output, networks, data_path, flood_areas, flood_layer):
    """Count assets from each network intersecting with flood polygons and save as a consolidated CSV."""
    output_dir = f"{output}/coastal_protection_assets"
    os.makedirs(output_dir, exist_ok=True)

    combined_results = []  # Accumulate all dataframes here

    for rcp in RCP:
        for rp in RP:
            # flood_layer = f"flood_protection_area_rcp_{rcp}_rp{rp}"
            flood_layer = flood_layer
            logging.info(f"Processing layer: {flood_layer}")

            try:
                # Load flood polygons
                flood_polygons = gpd.read_file(flood_areas, layer=flood_layer)

                if flood_polygons.empty or "id" not in flood_polygons.columns:
                    logging.warning(
                        f"No valid polygons or missing 'id' column in {flood_layer}"
                    )
                    continue

                # Initialize results dictionary for this scenario
                all_results = {str(pid): {} for pid in flood_polygons["id"]}

                cost_columns = []

                # Process each network
                for _, network in networks.iterrows():
                    ref = f"{network['asset_gpkg']}_{network['asset_layer']}"
                    cost_col = f"{ref}_cost"

                    network_results, cost_results = process_network_assets(
                        network, flood_polygons, data_path
                    )

                    for pid in all_results:
                        all_results[pid][ref] = network_results.get(pid, 0)
                        all_results[pid][cost_col] = round(
                            cost_results.get(pid, 0.0), 2
                        )

                # Convert results to DataFrame
                result_df = pd.DataFrame.from_dict(all_results, orient="index")
                result_df.index.name = "polygon_id"
                result_df.reset_index(inplace=True)

                # Ensure all expected asset columns exist
                unique_refs = [
                    f"{net['asset_gpkg']}_{net['asset_layer']}"
                    for _, net in networks.iterrows()
                ]
                for ref in unique_refs:
                    if ref not in result_df.columns:
                        result_df[ref] = 0
                    cost_col = f"{ref}_cost"
                    if cost_col not in result_df.columns:
                        result_df[cost_col] = 0.0

                cost_columns = [f"{ref}_cost" for ref in unique_refs]
                result_df["total_cost"] = result_df[cost_columns].sum(axis=1).round(2)

                # Extract epoch and rcp
                epoch = rcp[-4:]
                rcp_prefix = rcp[:-4]
                try:
                    rcp_value = float(rcp_prefix) / 10.0
                except ValueError:
                    rcp_value = rcp_prefix

                result_df["rcp"] = rcp_value
                result_df["epoch"] = epoch
                result_df["rp"] = rp

                combined_results.append(result_df)

                total_count = result_df[unique_refs].sum().sum()
                logging.info(
                    f"Processed layer {flood_layer} with {total_count} total intersections"
                )

            except Exception as e:
                logging.error(f"Error processing layer {flood_layer}: {str(e)}")

    if combined_results:
        final_df = pd.concat(combined_results, ignore_index=True)

        cost_columns = [col for col in final_df.columns if col.endswith("_cost")]
        other_columns = [col for col in final_df.columns if not col.endswith("_cost")]

        # Define your explicit order for non-cost columns
        explicit_order = ["polygon_id", "rcp", "epoch", "rp"]  # Your priority columns
        remaining_columns = [col for col in other_columns if col not in explicit_order]

        # Combine in desired order
        new_column_order = explicit_order + remaining_columns + cost_columns
        final_df = final_df[new_column_order]

        return final_df
    else:
        logging.warning("No data was processed, nothing to save.")


def add_coastline_lengths(RCP, RP, flood_asset_df, flood_areas, flood_layer):
    # Create a copy to avoid modifying the original dataframe
    updated_flood_asset_df = flood_asset_df.copy()

    # Initialize the coastline_length column if it doesn't exist
    if "coastline_length" not in updated_flood_asset_df.columns:
        updated_flood_asset_df["coastline_length"] = None

    # Convert merge columns to consistent types in the main dataframe
    updated_flood_asset_df["polygon_id"] = updated_flood_asset_df["polygon_id"].astype(
        str
    )
    updated_flood_asset_df["epoch"] = updated_flood_asset_df["epoch"].astype(str)
    updated_flood_asset_df["rcp"] = updated_flood_asset_df["rcp"].astype(str)
    updated_flood_asset_df["rp"] = updated_flood_asset_df["rp"].astype(str)

    for rcp in RCP:
        for rp in RP:
            # flood_layer = f"flood_protection_coastline_rcp_{rcp}_rp{rp}"
            flood_layer = flood_layer
            flood_polygons = gpd.read_file(flood_areas, layer=flood_layer)
            logging.info(f"Processing layer: {flood_layer}")

            epoch = rcp[-4:]
            rcp_prefix = rcp[:-4]
            try:
                rcp_value = float(rcp_prefix) / 10.0
            except ValueError:
                rcp_value = rcp_prefix

            # Add the matching columns to flood_polygons for merging
            flood_polygons["epoch"] = str(epoch)
            flood_polygons["rcp"] = str(rcp_value)
            flood_polygons["rp"] = str(rp)

            # Rename 'id' to 'polygon_id' and 'length' to 'coastline_length_temp'
            flood_polygons_renamed = flood_polygons.rename(
                columns={"id": "polygon_id", "length": "coastline_length_temp"}
            )

            # Convert all merge columns to string for consistent merging
            flood_polygons_renamed["polygon_id"] = flood_polygons_renamed[
                "polygon_id"
            ].astype(str)

            # Merge to find matches
            merged = updated_flood_asset_df.merge(
                flood_polygons_renamed[
                    ["polygon_id", "coastline_length_temp", "epoch", "rcp", "rp"]
                ],
                on=["polygon_id", "epoch", "rcp", "rp"],
                how="left",
            )

            # Update the coastline_length column where matches were found
            mask = merged["coastline_length_temp"].notna()
            updated_flood_asset_df.loc[mask, "coastline_length"] = merged.loc[
                mask, "coastline_length_temp"
            ].round()

            logging.info(f"Updated {mask.sum()} rows with coastline lengths")

    return updated_flood_asset_df


@click.command()
@click.version_option("1.0")
@click.option(
    "--network-csv",
    "-n",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to network layers csv",
)
@click.option(
    "--processed-data-path",
    "-d",
    required=True,
    type=click.Path(exists=False, dir_okay=True, file_okay=False, readable=True),
    help="Path to processed data",
)
@click.option(
    "--coastal-adaptation-assets",
    "-c",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to GPKG with coastal protection assets",
)
@click.option(
    "--rcp",
    "-rcp",
    required=True,
    type=float,
    help="RPS value",
)
@click.option(
    "--rp",
    "-rp",
    required=True,
    type=int,
    help="RP value",
)
@click.option(
    "--epoch",
    "-ep",
    required=True,
    type=int,
    help="epcoh value",
)
@click.option(
    "--output-dir",
    "-o",
    required=True,
    type=click.Path(exists=False, dir_okay=True, file_okay=False, readable=True),
    help="Path to output gpkg",
)
def main(
    network_csv,
    processed_data_path,
    coastal_adaptation_assets,
    rcp,
    epoch,
    rp,
    output_dir,
):

    RP = [f"{rp}"]

    fl_map_rcp = f"{int(rcp * 10)}{epoch}"
    RCP = [fl_map_rcp]

    data_path = processed_data_path
    flood_area_layer = "areas"
    flood_coastline_layer = "edges"

    networks_csv = network_csv
    networks = pd.read_csv(networks_csv)
    networks = networks[networks["asset_description"] != "buildings"]
    # print (networks)

    flood_areas = coastal_adaptation_assets
    output_file = f"{output_dir}/coastal_protection_assets/coastal_protection_assets_breakdown.csv"

    flood_asset_df = count_assets(
        RCP, RP, output_dir, networks, data_path, flood_areas, flood_area_layer
    )
    flood_asset_df = add_coastline_lengths(
        RCP, RP, flood_asset_df, flood_areas, flood_coastline_layer
    )

    # Save consolidated CSV
    flood_asset_df.to_csv(output_file, index=False)
    logging.info(
        f"Saved consolidated CSV with {len(flood_asset_df)} rows to {output_file}"
    )


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s",
        level=logging.INFO,
    )

    main()
