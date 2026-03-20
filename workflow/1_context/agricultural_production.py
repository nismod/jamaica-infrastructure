"""Assign agriculture GDP to land use layers in Jamaica"""

import os
import subprocess

import pandas as pd
import geopandas as gpd

from shapely.geometry import Point
import click
from jamaica_infrastructure.geo import (
    LOCAL_PROJ_CRS_EPSG,
    raster_rewrite,
    create_voronoi_layer,
    remove_geometry_collections,
)


@click.command()
@click.version_option("1.0")
@click.option(
    "--spam-path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, readable=True),
    help="Path to SPAM agriculture rasters directory.",
)
@click.option(
    "--spam-agriculture-outputs-path",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, writable=True),
    help="Path to SPAM agriculture outputs GeoPackage.",
)
@click.option(
    "--crop-details-path",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True),
    help="Path to crop details CSV.",
)
@click.option(
    "--land-use-path",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True),
    help="Path to Jamaica land use combined with sectors GeoPackage.",
)
@click.option(
    "--economic-output-path",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True),
    help="Path to economic output CSV.",
)
@click.option(
    "--output-areas",
    required=True,
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
    help="Path to output agriculture GDP GeoPackage.",
)
def main(
    spam_path,
    spam_agriculture_outputs_path,
    crop_details_path,
    land_use_path,
    economic_output_path,
    output_areas,
):
    """
    Process agricultural production data and allocate GDP to land use areas.

    Example usage:

        python workflow/1_context/agricultural_production.py \
            --spam-path incoming_data/agriculture_data \
            --spam-agriculture-outputs-path processed_data/agriculture_data/spam_agriculture_outputs.gpkg \
            --crop-details-path processed_data/agriculture_data/crop_details.csv \
            --land-use-path processed_data/land_type_and_use/jamaica_land_use_combined_with_sectors.gpkg \
            --economic-output-path processed_data/macroeconomic_data/NIP_2023.csv \
            --output-areas processed_data/agriculture_data/agriculture_gdp.gpkg

    This function performs the following steps:

    1. Reads and reprojects SPAM (Spatial Production Allocation Model) crop
       raster data for Jamaica, converting to CSV via gdal2xyz and merging into
       a single GeoDataFrame.
    2. Creates Voronoi polygons from crop point data for spatial area
       representation.
    3. Maps crop types to economic sector codes and computes production values
       in USD per square meter for each crop area.
    4. Reads agricultural land use data and performs a spatial join with crop
       yield areas to intersect land use polygons with crop data.
    5. Allocates daily GDP (converted from annual JMD millions) to agricultural
       land parcels, weighted by crop production value and area, using economic
       output data filtered by sector code 'A'.
    6. Handles special cases: Post-harvest output (subsector 14) allocated by
       crop tonnage. Forest output (subsector 20) allocated by forest area.
       Non-agricultural land use areas are assigned zero GDP.
    7. Computes per-square-meter GDP and crop tonnage densities for each land
       parcel.
    9. Writes the following outputs: SPAM agriculture outputs - production_value,
       production_areas, tonnage_areas, yield_areas, values_per_sqm with crop point
       and area data; and areas GeoPackage with GDP and crop tonnage by land use
       area.
    """

    crop_folders = [
        "spam2010v2r0_global_val_prod_agg.geotiff",
        "spam2010v2r0_global_prod.geotiff",
        "spam2010v2r0_global_yield.geotiff",
    ]
    crop_strings = [
        "spam2010V2r0_global_V_agg_",
        "spam2010V2r0_global_P_",
        "spam2010V2r0_global_Y_",
    ]
    crop_outputs = ["production", "tonnage", "yield"]
    for crop_path, crop_field, crop_layer in zip(
        crop_folders, crop_strings, crop_outputs
    ):
        crop_data_path = os.path.join(spam_path, crop_path, "JAM")
        all_crops = []
        all_fields = []
        for file in os.listdir(crop_data_path):
            if file.endswith(".tif"):
                crop_file = file.replace(".tif", "")
                field_name = crop_file.replace(crop_field, "")
                all_fields.append(field_name)
                crop_raster_in = os.path.join(crop_data_path, f"{crop_file}.tif")
                crop_raster_out = os.path.join(
                    crop_data_path, f"{crop_file}_reproject.tif"
                )
                outCSVName = os.path.join(crop_data_path, f"{crop_file}.csv")
                if not os.path.exists(outCSVName):
                    raster_rewrite(crop_raster_in, crop_raster_out)
                    subprocess.run(["gdal2xyz.py", "-csv", crop_raster_out, outCSVName])

                # Load points and convert to geodataframe with coordinates
                load_points = pd.read_csv(
                    outCSVName,
                    header=None,
                    names=["x", "y", field_name],
                    index_col=None,
                )

                if len(all_crops) > 0:
                    all_crops = pd.merge(
                        all_crops, load_points, how="left", on=["x", "y"]
                    )
                else:
                    all_crops = load_points.copy()

                del load_points
                print("* Done with", file)

        all_crops["geometry"] = [Point(xy) for xy in zip(all_crops.x, all_crops.y)]
        crop_points = gpd.GeoDataFrame(
            all_crops, crs=f"EPSG:{LOCAL_PROJ_CRS_EPSG}", geometry="geometry"
        )
        crop_points["crop_id"] = crop_points.index.values.tolist()
        del all_crops

        crop_areas = create_voronoi_layer(
            crop_points, "crop_id", epsg=LOCAL_PROJ_CRS_EPSG
        )

        crop_areas = gpd.GeoDataFrame(
            pd.merge(
                crop_areas,
                crop_points[["crop_id"] + all_fields],
                how="left",
                on=["crop_id"],
            ),
            geometry="geometry",
            crs=f"EPSG:{LOCAL_PROJ_CRS_EPSG}",
        )

        crop_points.to_file(
            spam_agriculture_outputs_path,
            layer=f"{crop_layer}_value",
            driver="GPKG",
        )
        crop_areas.to_file(
            spam_agriculture_outputs_path,
            layer=f"{crop_layer}_areas",
            driver="GPKG",
        )
        del crop_areas, crop_points

    crop_yields = gpd.read_file(
        spam_agriculture_outputs_path,
        layer=f"tonnage_areas",
    )
    crop_yields = crop_yields.to_crs(epsg=LOCAL_PROJ_CRS_EPSG)
    crop_details = pd.read_csv(crop_details_path)
    tech_type = ["A", "I", "R"]
    poultry_crops = ["maiz", "ocer", "pmil", "smil", "soyb", "sunf", "whea"]

    crop_yields["crop_tons"] = 0

    all_crop_columns = []
    all_sector_columns = []
    for crop in crop_details.itertuples():
        crop_columns = [f"{crop.name.upper()}_{t}" for t in tech_type]
        sector_columns = [
            f"{crop.sector_code}_{crop.subsector_code}_{t}" for t in tech_type
        ]
        poultry_columns = [f"A_12_{t}" for t in tech_type]
        all_sector_columns += sector_columns + poultry_columns
        for i, (cr, sc, pc) in enumerate(
            list(zip(crop_columns, sector_columns, poultry_columns))
        ):
            if sc not in crop_yields.columns.values.tolist():
                crop_yields[sc] = 0
            if pc not in crop_yields.columns.values.tolist():
                crop_yields[pc] = 0
            crop_yields[cr] = crop_yields.apply(
                lambda x: x[cr] if x[cr] > 0 else 0, axis=1
            )
            crop_yields[f"{cr}_prod"] = crop.value_usd_per_ton * crop_yields[cr]
            crop_yields["crop_tons"] += crop_yields[cr]
            crop_yields[sc] += crop_yields[f"{cr}_prod"]
            if crop.name in poultry_crops:
                crop_yields[pc] += crop_yields[f"{cr}_prod"]
            all_crop_columns.append(f"{cr}_prod")

        print("* Done with", crop.name, crop.value_usd_per_ton)

    crop_yields["crop_tons_persqm"] = crop_yields["crop_tons"] / crop_yields["areas"]
    crop_yields = crop_yields[
        ["crop_id", "areas", "crop_tons_persqm", "geometry"]
        + list(set(all_sector_columns))
        + all_crop_columns
    ]

    crop_yields.to_file(
        spam_agriculture_outputs_path,
        layer=f"values_per_sqm",
        driver="GPKG",
    )

    #
    # Read economic data
    #
    economic_output = pd.read_csv(economic_output_path, comment="#").query(
        "Code == 'A'"
    )
    value_colname = "GVA JMD Millions (2023)"

    agri_land_use = gpd.read_file(
        land_use_path,
        layer="areas",
    )
    agri_land_use = agri_land_use.to_crs(epsg=LOCAL_PROJ_CRS_EPSG)
    agri_land_use["land_id"] = agri_land_use.index.values.tolist()
    agri_land_use["land_id"] = agri_land_use.progress_apply(
        lambda x: f"land_{x.land_id}", axis=1
    )
    non_agri_land_use = agri_land_use[
        ~(
            (agri_land_use["sector_code_forest"] == "A")
            | (agri_land_use["sector_code_tnc"] == "A")
        )
    ]
    non_agri_land_use["A_GDP"] = 0
    agri_land_use = agri_land_use[
        (agri_land_use["sector_code_forest"] == "A")
        | (agri_land_use["sector_code_tnc"] == "A")
    ]
    # TODO handle subsectors
    agri_land_use["known_forest"] = agri_land_use.progress_apply(
        lambda x: (
            1
            if (
                str(x.subsector_code_forest) == "20"
                and str(x.subsector_code_tnc) == "20"
            )
            else 0
        ),
        axis=1,
    )
    agri_forest = agri_land_use[agri_land_use["known_forest"] == 1]
    agri_land_use = agri_land_use[agri_land_use["known_forest"] == 0]
    for del_col in ["index", "index_left", "index_right"]:
        if del_col in agri_land_use.columns.values.tolist():
            agri_land_use.drop(del_col, axis=1, inplace=True)

    agri_areas = gpd.sjoin(
        agri_land_use, crop_yields, how="inner", predicate="intersects"
    ).reset_index()
    agri_areas.rename(columns={"geometry": "agri_geometry"}, inplace=True)
    agri_areas = pd.merge(
        agri_areas, crop_yields[["crop_id", "geometry"]], how="left", on=["crop_id"]
    )
    agri_areas["geom"] = agri_areas.progress_apply(
        lambda x: x["agri_geometry"].intersection(x["geometry"].buffer(0)), axis=1
    )
    agri_areas.drop(["agri_geometry", "geometry"], axis=1, inplace=True)
    agri_areas.rename(columns={"geom": "geometry"}, inplace=True)
    values_columns = (
        [
            "land_id",
            "forest_id",
            "Classify",
            "LU_CODE",
            "tnc_id",
            "NAME",
            "TNCCODE",
            "global_id",
            "global_LU_type",
            "sector_code_forest",
            "sector_code_tnc",
            "crop_tons_persqm",
        ]
        + list(set(all_sector_columns))
        + all_crop_columns
        + ["geometry"]
    )
    agri_areas = agri_areas[values_columns]
    agri_areas["area_sqm"] = agri_areas.progress_apply(
        lambda x: x.geometry.area, axis=1
    )
    agri_areas["A_GDP"] = 0
    tot_gpd = 0
    for i, econ in economic_output.iterrows():
        econ_subsector_codes = str(econ["subsector_code"]).split(",")
        econ_codes = [
            f"A_{e}_A"
            for e in econ_subsector_codes
            if f"A_{e}_A" in agri_areas.columns.values.tolist()
        ]
        if len(econ_codes) > 0:
            output = (1.0e6 / 365.0) * econ[value_colname]
            tot_gpd += output
            agri_areas["weight"] = (
                agri_areas[econ_codes].sum(axis=1) * agri_areas["area_sqm"]
            )
            agri_areas["A_GDP"] += (
                output * agri_areas["weight"] / agri_areas["weight"].sum()
            )

    post_harvest_output = (1.0e6 / 365.0) * economic_output[
        (economic_output["Code"] == "A") & (economic_output["subsector_code"] == "14")
    ][value_colname].sum()
    tot_gpd += post_harvest_output
    print("Sector A-14 output", post_harvest_output)

    agri_areas["crop_tons"] = agri_areas["crop_tons_persqm"] * agri_areas["area_sqm"]
    agri_areas["A_GDP"] += (
        post_harvest_output * agri_areas["crop_tons"] / agri_areas["crop_tons"].sum()
    )

    agri_areas = (
        agri_areas.groupby(["land_id"])["A_GDP", "crop_tons"].sum().reset_index()
    )
    agri_areas = pd.merge(agri_land_use, agri_areas, how="left", on=["land_id"])

    forest_output = (1.0e6 / 365.0) * economic_output[
        (economic_output["Code"] == "A") & (economic_output["subsector_code"] == "20")
    ][value_colname].sum()
    tot_gpd += forest_output
    agri_forest["A_GDP"] = (
        forest_output * agri_forest["area_m2"] / agri_forest["area_m2"].sum()
    )
    values_columns = [
        "land_id",
        "forest_id",
        "Classify",
        "LU_CODE",
        "tnc_id",
        "NAME",
        "TNCCODE",
        "global_id",
        "global_LU_type",
        "sector_code_forest",
        "sector_code_tnc",
        "area_m2",
        "A_GDP",
        "geometry",
    ]
    agri_areas = pd.concat(
        [
            agri_areas[values_columns + ["crop_tons"]],
            agri_forest[values_columns],
            non_agri_land_use[values_columns],
        ],
        axis=0,
        ignore_index=True,
    )
    agri_areas["crop_tons"] = agri_areas["crop_tons"].fillna(0)
    agri_areas["A_GDP"] = agri_areas["A_GDP"].fillna(0)

    agri_areas["GDP_persqm"] = agri_areas["A_GDP"] / agri_areas["area_m2"]
    agri_areas["crop_tons_persqm"] = agri_areas["crop_tons"] / agri_areas["area_m2"]
    agri_areas["GDP_unit"] = "JD/day"

    agri_areas = gpd.GeoDataFrame(
        agri_areas, geometry="geometry", crs=f"EPSG:{LOCAL_PROJ_CRS_EPSG}"
    )
    agri_areas = remove_geometry_collections(agri_areas)

    tot_area = agri_areas["area_m2"].sum()
    print("* Given GDP", tot_gpd)
    print("* Estimated GDP", agri_areas["A_GDP"].sum())
    print("* Estimated Areas", tot_area)

    gpd.GeoDataFrame(
        agri_areas, geometry="geometry", crs=f"EPSG:{LOCAL_PROJ_CRS_EPSG}"
    ).to_file(
        output_areas,
        layer="areas",
        driver="GPKG",
    )


if __name__ == "__main__":
    # TODO logging.basicConfig
    main()
