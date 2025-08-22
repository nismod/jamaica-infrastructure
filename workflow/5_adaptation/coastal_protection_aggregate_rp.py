"""
Calculate the per return period damages and losses at a coastal defence feature
level (i.e. total across assets for each sea wall segment)
"""


import logging
from pathlib import Path

import click
import numpy as np
import pandas as pd

from jamaica_infrastructure.utils import is_sole_value


def read_lookup(path: str) -> pd.DataFrame:
    return pd.read_parquet(path) \
        .rename(columns={"flood_id_rcp_852100_rp_100": "defence_id"}) \
        .astype({"defence_id": int})


def read_risk(
    asset_id_column: str,
    return_period_risk_path: str,
    kind: str,
    lookup_path: str,
) -> pd.DataFrame:

    per_layer_risk = pd.read_parquet(return_period_risk_path).set_index(asset_id_column)

    # We will drop and then rewrite the units later; ensure they're as expected
    if kind == "damage":
        assert set(per_layer_risk["damage_cost_unit"].unique()).issubset({"J$", "JD"})
    elif kind == "loss":
        assert is_sole_value(per_layer_risk["economic_loss_unit"], "J$/day")
    else:
        raise ValueError(f"{kind=}, but must be either 'damage' or 'loss'")

    risk_cols = [c for c in per_layer_risk.columns if c.startswith("coastal__rp_")]
    per_layer_risk = per_layer_risk.loc[:, risk_cols]

    try:
        lookup = read_lookup(lookup_path)
    except TypeError:
        logging.info(f".. Failed to read {lookup_path=}")
        return pd.DataFrame([])

    return lookup \
        .set_index(asset_id_column) \
        .loc[:, ["defence_id"]] \
        .join(per_layer_risk, how="left") \
        .reset_index() \
        .rename(columns={asset_id_column: "asset_id"})


@click.command()
@click.version_option("1.0")
@click.option(
    "--damage",
    "-d",
    "damage_paths",
    required=True,
    multiple=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Damage (per return period) files (one per layer)",
)
@click.option(
    "--loss",
    "-l",
    "loss_paths",
    required=True,
    multiple=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Economic loss (per return period) files (one per layer)",
)
@click.option(
    "--lookup",
    "-lk",
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
    "--output-exposure",
    "-oe",
    "output_exposure_path",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True, readable=True),
    help="Location to write per coastal protection feature exposure",
)
@click.option(
    "--output-damage",
    "-od",
    "output_damage_path",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True, readable=True),
    help="Location to write per coastal protection feature avoided damages",
)
@click.option(
    "--output-loss",
    "-ol",
    "output_loss_path",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True, readable=True),
    help="Location to write per coastal protection feature avoided losses",
)
def coastal_protection_aggregate(
    damage_paths: tuple[str],
    loss_paths: tuple[str],
    lookup_paths: tuple[str],
    network_path: str,
    output_exposure_path: str,
    output_damage_path: str,
    output_loss_path: str,
) -> None:

    network: pd.DataFrame = pd.read_csv(network_path)

    damages_by_layer = []
    losses_by_layer = []
    for row in network.itertuples():

        gpkg_layer: str = f"{row.asset_gpkg}_{row.asset_layer}"
        logging.info(f"{gpkg_layer}")

        damage_path, = filter(lambda x: Path(x).name == f"{gpkg_layer}_damages.parquet", damage_paths)
        loss_path, = filter(lambda x: Path(x).name == f"{gpkg_layer}_losses.parquet", loss_paths)
        lookup_path, = filter(lambda x: Path(x).name == f"{gpkg_layer}_coastal_filtered.parquet", lookup_paths)

        damages_by_layer.append(read_risk(row.asset_id_column, damage_path, "damage", lookup_path))
        losses_by_layer.append(read_risk(row.asset_id_column, loss_path, "loss", lookup_path))

    damages = pd.concat(damages_by_layer).rename(columns={"defence_id": "id"})
    damages = damages.drop(columns=["asset_id"]).groupby("id").sum().reset_index()
    damages["damage_cost_unit"] = "J$"
    damages.to_parquet(output_damage_path)

    losses = pd.concat(losses_by_layer).rename(columns={"defence_id": "id"})
    losses = losses.drop(columns=["asset_id"]).groupby("id").sum().reset_index()
    exposures = losses.set_index("id").copy()
    losses["economic_loss_unit"] = "J$/day"
    losses.to_parquet(output_loss_path)

    # We can't produce an aggregate exposure for each coastal protection feature,
    # as the protected assets use a variety of units (no unit, m, m2, etc.)
    # We write an empty file to conform with the existing ETL process
    exposures.loc[:, :] = np.nan
    exposures["exposure_unit"] = ""
    exposures.reset_index().to_parquet(output_exposure_path)


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    coastal_protection_aggregate()
