"""
Associate power system failure with economic losses.

Reads disruption results for all grid cells and joins with GDP data
at affected nodes to estimate economic losses.
"""

import logging
from pathlib import Path

import click
import pandas as pd
from tqdm.auto import tqdm


@click.command()
@click.version_option("1.0.0")
@click.option(
    "--node-gdp-path",
    "-g",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Path to CSV with node GDP data.",
)
@click.option(
    "--disruption-dir",
    "-d",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True),
    help="Directory containing disruption CSV files.",
)
@click.option(
    "--output-dir",
    "-o",
    required=True,
    type=click.Path(file_okay=False, dir_okay=True, writable=True),
    help="Output directory for loss CSV files.",
)
def main(node_gdp_path, disruption_dir, output_dir):
    """Calculate economic losses from grid cell disruptions."""
    logging.info(f"Reading node GDP data from {node_gdp_path}")
    node_gdp = pd.read_csv(
        node_gdp_path,
        index_col="id",
        usecols=["id", "total_GDP"],
    )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    disruption_dir = Path(disruption_dir)
    csv_files = list(disruption_dir.glob("disruption_*.csv"))
    logging.info(f"Processing {len(csv_files)} disruption files")

    for csv_path in tqdm(csv_files):
        disruption = pd.read_csv(csv_path)
        disruption.set_index("affected_node_id", inplace=True)

        # Get cell coordinates from the disruption data
        (cell_x,) = disruption.cell_index_x.unique()
        (cell_y,) = disruption.cell_index_y.unique()

        disruption_loss = disruption.join(node_gdp).rename(
            columns={"total_GDP": "loss_gdp"}
        )
        disruption_loss["loss_gdp_unit"] = "JD/day"

        output_path = output_dir / f"loss_{cell_x}_{cell_y}.csv"
        disruption_loss.to_csv(output_path)

    logging.info(f"Wrote {len(csv_files)} loss files to {output_dir}")


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    main()
