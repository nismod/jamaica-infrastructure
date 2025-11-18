"""
Calculate the Benefit Cost Ratio (BCR) per coastal defence feature (e.g. sea
wall segment) where the cost is the cost of the segment and the benefit is the
Net Present Value (NPV) of the avoided damage and loss.
"""

import logging
from pathlib import Path
import re

import click
import pandas as pd


def read_lookup(path: str) -> pd.DataFrame:
    return (
        pd.read_parquet(path).rename(columns={"flood_id_rcp_852100_rp_100": "defence_id"}).astype({"defence_id": int})
    )


def collate_costs_and_avoided_risks(
    gpkg_layer: str,
    asset_id_column: str,
    EAD_EAEL_path: str,
    adaptation_cost_path: str,
    lookup_path: str,
) -> pd.DataFrame:

    per_layer_EAD_EAEL = pd.read_csv(EAD_EAEL_path).drop(columns=["adapt_cost_npv"]).set_index(asset_id_column)

    per_layer_costs = pd.read_csv(adaptation_cost_path)
    if per_layer_costs.empty:
        logging.info(f".. No cost data for {gpkg_layer}")
        return pd.DataFrame([])

    try:
        lookup = read_lookup(lookup_path)
    except TypeError:
        logging.info(f".. Failed to read lookup for {gpkg_layer}")
        return pd.DataFrame([])

    return (
        lookup.set_index(asset_id_column)
        .loc[:, ["defence_id"]]
        .join(per_layer_EAD_EAEL, how="left")
        .join(per_layer_costs.set_index(asset_id_column).loc[:, "adapt_cost_npv"])
        .reset_index()
        .rename(columns={asset_id_column: "asset_id"})
    )


@click.command()
@click.version_option("1.0")
@click.option(
    "--cost",
    "-c",
    "adaptation_cost_paths",
    required=True,
    multiple=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Adaptation option cost files (per layer)",
)
@click.option(
    "--EAD-EAEL",
    "-ee",
    "EAD_EAEL_paths",
    required=True,
    multiple=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="EAD and EAEL files (per layer)",
)
@click.option(
    "--lookup",
    "-l",
    "lookup_paths",
    required=True,
    multiple=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Lookup table from coastal defence asset to protected assets (per layer)",
)
@click.option(
    "--network",
    "-n",
    "network_path",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True, readable=True),
    help="Path to network layers table",
)
@click.option(
    "--output-BCR",
    "-ob",
    "output_BCR_path",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True, readable=True),
    help="Location to write per coastal protection asset BCR",
)
@click.option(
    "--output-EAD-EAEL",
    "-oe",
    "output_EAD_EAEL_path",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True, readable=True),
    help="Location to write per coastal protection asset avoided EAD & EAEL",
)
def coastal_protection_aggregate(
    adaptation_cost_paths: tuple[str],
    EAD_EAEL_paths: tuple[str],
    lookup_paths: tuple[str],
    network_path: str,
    output_BCR_path: str,
    output_EAD_EAEL_path: str,
) -> None:

    network: pd.DataFrame = pd.read_csv(network_path)

    by_layer = []
    for row in network.itertuples():

        gpkg_layer: str = f"{row.asset_gpkg}_{row.asset_layer}"
        logging.info(f"{gpkg_layer}")

        (EAD_EAEL_path,) = filter(
            lambda x: Path(x).name == f"coastal_{gpkg_layer}_adaptation_costs_avoided_EAD_EAEL.csv", EAD_EAEL_paths
        )
        (adaptation_cost_path,) = filter(
            lambda x: Path(x).name == f"{gpkg_layer}_adaptation_timeseries_and_npvs.csv", adaptation_cost_paths
        )
        (lookup_path,) = filter(lambda x: Path(x).name == f"{gpkg_layer}_coastal_filtered.parquet", lookup_paths)

        by_layer.append(
            collate_costs_and_avoided_risks(
                gpkg_layer, row.asset_id_column, EAD_EAEL_path, adaptation_cost_path, lookup_path
            )
        )

    data = pd.concat(by_layer).rename(columns={"defence_id": "id"})
    (adaptation_option_name,) = data["adaptation_option"].dropna().unique()

    old_risk_cols = []
    rcps = set()
    for col in data.columns:
        if match := re.match(r"coastal__rcp_(\d\.\d)__avoided_(EAD|EAEL)_(amin|mean|amax)", col):
            rcp, var, agg = match.groups()
            old_risk_cols.append(col)
            rcps.add(rcp)

    # Frontend assumes hazard={'TC', 'flooding'}, so we rename here
    new_risk_cols = [c.replace("coastal", "flooding") for c in old_risk_cols]
    data = data.rename(columns=dict(zip(old_risk_cols, new_risk_cols)))

    per_segment_costs = data.loc[:, ["id", "adapt_cost_npv"]].drop_duplicates().copy()
    per_segment_costs = (
        per_segment_costs[~per_segment_costs["adapt_cost_npv"].isna()].sort_values("id").reset_index(drop=True)
    )
    assert len(per_segment_costs) == len(per_segment_costs["id"].unique())

    # Find the sum (total) avoided risk per coastal defence segment
    EAD_EAEL = data.loc[:, ["id"] + new_risk_cols].groupby("id").sum()
    EAD_EAEL = EAD_EAEL.join(per_segment_costs.set_index("id"))
    EAD_EAEL["adaptation_option"] = adaptation_option_name
    EAD_EAEL["adapt_cost_unit"] = "J$"
    EAD_EAEL["ead_cost_unit"] = "J$"
    EAD_EAEL["eael_cost_unit"] = "J$/day"
    # We used an n year return period map to define flood defence height
    EAD_EAEL["flood_depth_protection_level"] = 100
    EAD_EAEL.to_csv(output_EAD_EAEL_path)

    # Calculate the Benefit Cost Ratio (total avoided risk / adaptation cost)
    BCR = EAD_EAEL.loc[:, ["adapt_cost_npv"] + new_risk_cols].copy()
    for rcp in sorted(rcps):
        for var in ("amin", "mean", "amax"):
            BCR[f"flooding__rcp_{rcp}__avoided_risk_{var}"] = (
                BCR[f"flooding__rcp_{rcp}__avoided_EAD_{var}"] + BCR[f"flooding__rcp_{rcp}__avoided_EAEL_{var}"]
            )
            BCR[f"flooding__rcp_{rcp}__BCR_{var}"] = (
                BCR[f"flooding__rcp_{rcp}__avoided_risk_{var}"] / BCR["adapt_cost_npv"]
            )
    BCR = BCR.drop(columns=new_risk_cols)
    BCR["adaptation_option"] = adaptation_option_name
    BCR["cost_units"] = "J$"
    BCR["flood_depth_protection_level"] = 100
    BCR.to_csv(output_BCR_path)

    return


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    coastal_protection_aggregate()
