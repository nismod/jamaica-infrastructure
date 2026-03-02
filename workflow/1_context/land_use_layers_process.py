"""Create a combined land use layer from TNC and Forestry land use layers
    Also combine a global mining areas land use layer 
    Add macroeconomic sector and subsector codes to the data
"""

import sys
import os

import pandas as pd
import geopandas as gpd
import numpy as np
from preprocess_utils import *
from tqdm import tqdm
import click

from jamaica_infrastructure import LOCAL_PROJ_CRS_EPSG

tqdm.pandas()

@click.command()
@click.version_option("1.0")
@click.option(
    "--data-dir",
    "-d",
    required=True,
    type=click.Path(exists=False, dir_okay=True, file_okay=False, readable=True),
    help="Path to processed data",
)
@click.option(
    "--incoming-data-dir",
    "-i",
    required=True,
    type=click.Path(exists=False, dir_okay=True, file_okay=False, readable=True),
    help="Path to unprocessed incoming data",
)
@click.option(
    "--tnc-landuse-path",
    "-tl",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="",
)
@click.option(
    "--forest-landuse-path",
    "-fl",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="",
)
@click.option(
    "--mining-landuse-path",
    "-ml",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="",
)
@click.option(
    "--forest-sector-mapping-path",
    "-fsm",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="",
)
@click.option(
    "--tnc-sector-mapping-path",
    "-tsm",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="",
)

def main(
    data_dir,
    incoming_data_dir,
    tnc_landuse_path,
    forest_landuse_path,
    mining_landuse_path,
    forest_sector_mapping_path,
    tnc_sector_mapping_path,
):
    incoming_data_path = incoming_data_dir
    processed_data_path = data_dir
    # database_name = "GWP_Jamaica_NSP_Master_Geodatabase_v01.gdb"

    tnc_landuse = gpd.read_file(
        tnc_landuse_path,
        layer="LandUse_LandUse",
    ).to_crs(epsg=LOCAL_PROJ_CRS_EPSG)
    tnc_landuse["tnc_id"] = tnc_landuse.index.values.tolist()
    forest_landuse = gpd.read_file(
        forest_landuse_path
    ).to_crs(epsg=LOCAL_PROJ_CRS_EPSG)
    forest_landuse["forest_id"] = forest_landuse.index.values.tolist()
    mining_landuse = gpd.read_file(
        mining_landuse_path
    ).to_crs(epsg=LOCAL_PROJ_CRS_EPSG)
    mining_landuse = mining_landuse[mining_landuse["COUNTRY_NAME"] == "Jamaica"]
    mining_landuse["global_id"] = mining_landuse.index.values.tolist()
    mining_landuse["global_LU_type"] = "Bauxite Extraction"

    tnc_landuse.to_file(
        os.path.join(
            processed_data_path, "land_type_and_use", "input_land_use_layers.gpkg"
        ),
        layer="tnc_land_use",
        driver="GPKG",
    )
    forest_landuse.to_file(
        os.path.join(
            processed_data_path, "land_type_and_use", "input_land_use_layers.gpkg"
        ),
        layer="forest_land_use",
        driver="GPKG",
    )
    mining_landuse.to_file(
        os.path.join(
            processed_data_path, "land_type_and_use", "input_land_use_layers.gpkg"
        ),
        layer="global_mining_land_use",
        driver="GPKG",
    )
    """
    Intersect the global mining and TNC land use layers and merge the modified geometries
    """
    tnc_global_mining_df = spatial_scenario_selection(
        tnc_landuse,
        mining_landuse,
        ["tnc_id", "NAME", "TNCCODE"],
        ["global_id", "global_LU_type"],
    )
    tnc_global_mining_df = gpd.GeoDataFrame(
        pd.DataFrame(tnc_global_mining_df),
        geometry="geometry",
        crs=f"EPSG:{LOCAL_PROJ_CRS_EPSG}",
    )
    merge_geometry = tnc_global_mining_df.dissolve(by="tnc_id").reset_index()
    tnc_r = []
    for values in merge_geometry.itertuples():
        val_id = getattr(values, "tnc_id")
        val_name = getattr(values, "NAME")
        val_code = getattr(values, "TNCCODE")
        tnc_remaining = (
            tnc_landuse[tnc_landuse["tnc_id"] == val_id]["geometry"]
            .values[0]
            .difference(values.geometry)
        )
        if tnc_remaining.is_empty is False:
            tnc_r.append((val_id, val_name, val_code, tnc_remaining))

    tnc_r = pd.DataFrame(tnc_r, columns=["tnc_id", "NAME", "TNCCODE", "geometry"])
    tnc_modified = tnc_landuse[
        ~tnc_landuse["tnc_id"].isin(tnc_r.tnc_id.values.tolist())
    ][["tnc_id", "NAME", "TNCCODE", "geometry"]]
    tnc_final = gpd.GeoDataFrame(
        pd.concat(
            [tnc_global_mining_df, tnc_r, tnc_modified], axis=0, ignore_index=True
        ),
        geometry="geometry",
        crs=f"EPSG:{LOCAL_PROJ_CRS_EPSG}",
    )
    tnc_final.to_file(
        os.path.join(
            processed_data_path, "land_type_and_use", "land_use_modified.gpkg"
        ),
        layer="tnc_modified",
        driver="GPKG",
    )

    """
    Intersect the modified TNC layer and the forest layer to create the final layer and geometries
    """
    tnc_final["tnc_geomtery"] = tnc_final.geometry
    forest_landuse["forest_geomtery"] = forest_landuse.geometry
    matches = gpd.sjoin(
        forest_landuse[
            ["forest_id", "Classify", "LU_CODE", "forest_geomtery", "geometry"]
        ],
        tnc_final[
            [
                "tnc_id",
                "NAME",
                "TNCCODE",
                "global_id",
                "global_LU_type",
                "tnc_geomtery",
                "geometry",
            ]
        ],
        how="inner",
        op="intersects",
    ).reset_index()
    matches = matches.drop("geometry", axis=1)
    matches["geometry"] = matches.progress_apply(
        lambda x: (x.tnc_geomtery.buffer(0)).intersection(x.forest_geomtery.buffer(0)),
        axis=1,
    )
    matches = matches.drop(["tnc_geomtery", "forest_geomtery"], axis=1)
    tnc_forest_df = gpd.GeoDataFrame(
        matches[~matches.is_empty], geometry="geometry", crs=f"EPSG:{LOCAL_PROJ_CRS_EPSG}"
    )

    """If we want to save the intermediate result
    """
    # tnc_forest_df.to_file(
    #                     os.path.join(
    #                     processed_data_path,
    #                     "land_type_and_use",
    #                     "land_use_combined.gpkg"
    #                         ),
    #                     layer="land_use_modified",
    #                     driver="GPKG"
    #                 )
    # tnc_forest_df = gpd.read_file(os.path.join(
    #                     processed_data_path,
    #                     "land_type_and_use",
    #                     "land_use_combined.gpkg"
    #                         ),
    #                     layer="land_use_modified")

    tnc_forest_df = split_multigeometry(tnc_forest_df)
    tnc_forest_df = gpd.GeoDataFrame(
        tnc_forest_df[
            ~(tnc_forest_df.geometry.geom_type.isin(["LineString", "Point"]))
        ],
        geometry="geometry",
        crs=f"EPSG:{LOCAL_PROJ_CRS_EPSG}",
    )
    tnc_forest_df["area_m2"] = tnc_forest_df.progress_apply(
        lambda x: x.geometry.area, axis=1
    )
    tnc_forest_df["area_hectare"] = tnc_forest_df["area_m2"] / 10000.0
    tnc_forest_df.to_file(
        os.path.join(
            processed_data_path, "land_type_and_use", "jamaica_land_use_combined.gpkg"
        ),
        layer="areas",
        driver="GPKG",
    )
    tnc_forest_df.groupby("Classify")["area_m2"].sum().reset_index().to_csv(
        os.path.join(processed_data_path, "land_type_and_use", "forest_classes.csv"),
        index=False,
    )
    tnc_forest_df.groupby("NAME")["area_m2"].sum().reset_index().to_csv(
        os.path.join(processed_data_path, "land_type_and_use", "tnc_classes.csv"),
        index=False,
    )

    """Add the economics sector mapping information
    """
    tnc_forest_df = gpd.read_file(
        os.path.join(
            processed_data_path, "land_type_and_use", "jamaica_land_use_combined.gpkg"
        ),
        layer="areas",
    )
    forest_sector_mapping = pd.read_csv(
        forest_sector_mapping_path
    )
    forest_sector_mapping.rename(
        columns={
            "sector_code": "sector_code_forest",
            "subsector_code": "subsector_code_forest",
            "infra_sector": "infra_sector_forest",
        },
        inplace=True,
    )
    tnc_sector_mapping = pd.read_csv(
        tnc_sector_mapping_path
    )
    tnc_sector_mapping.rename(
        columns={
            "sector_code": "sector_code_tnc",
            "subsector_code": "subsector_code_tnc",
            "infra_sector": "infra_sector_tnc",
        },
        inplace=True,
    )

    tnc_forest_df = pd.merge(
        tnc_forest_df, forest_sector_mapping, how="left", on=["Classify"]
    )
    tnc_forest_df = pd.merge(tnc_forest_df, tnc_sector_mapping, how="left", on=["NAME"])

    gpd.GeoDataFrame(
        tnc_forest_df, geometry="geometry", crs=f"EPSG:{LOCAL_PROJ_CRS_EPSG}"
    ).to_file(
        os.path.join(
            processed_data_path,
            "land_type_and_use",
            "jamaica_land_use_combined_with_sectors.gpkg",
        ),
        layer="areas",
        driver="GPKG",
    )


if __name__ == "__main__":
    main()
