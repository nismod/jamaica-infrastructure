"""Calculate proportions of other land cover classes within a buffer strip
- around various groups of land cover classes
- for various buffer distances

Example usage:

python surrounding_landcover.py \
    "~/data/jamaica-local/incoming/Terrestrial Land Cover (From Forestry Department)" \
    ~/data/jamaica-local/results/nbs \
    secondary-forest \
    100
"""

import os
import sys
import warnings

import geopandas
import pandas
import shapely
from tqdm import tqdm


CLASS_GROUPS = {
    "primary-forest": ["Closed broadleaved forest (Primary Forest)"],
    "secondary-forest": [
        "Bamboo and Secondary Forest",
        "Disturbed broadleaved forest (Secondary Forest)",
        "Fields and Secondary Forest",
        "Fields or Secondary Forest/Pine Plantation",
        "Secondary Forest",
    ],
    "dry-forest": [
        "Open dry forest - Short",
        "Open dry forest - Tall (Woodland/Savanna)",
    ],
    "wetland": ["Herbaceous Wetland"],
    "plantations": [
        "Plantation: Tree crops, shrub crops, sugar cane, banana",
        "Hardwood Plantation: Euculytus",
        "Hardwood Plantation: Mahoe",
        "Hardwood Plantation: Mahogany",
        "Hardwood Plantation: Mixed",
    ],
    "fields": [
        "Bamboo and Fields",
        "Fields  and Bamboo",
        "Fields: Bare Land",
        "Fields: Herbaceous crops, fallow, cultivated vegetables",
        "Fields: Pasture,Human disturbed, grassland",
    ],
    "agriculture": [
        "Bamboo and Fields",
        "Fields  and Bamboo",
        "Fields: Bare Land",
        "Fields: Herbaceous crops, fallow, cultivated vegetables",
        "Fields: Pasture,Human disturbed, grassland" "Hardwood Plantation: Euculytus",
        "Hardwood Plantation: Mahoe",
        "Hardwood Plantation: Mahogany",
        "Hardwood Plantation: Mixed"
        "Plantation: Tree crops, shrub crops, sugar cane, banana",
    ],
}


def main(input_path, output_path, class_group, buffer_distance):
    output_csv = os.path.join(
        output_path,
        f"2013_landuse_landcover__{class_group}__within_{buffer_distance}m.csv",
    )
    output_gpkg = os.path.join(
        output_path,
        f"2013_landuse_landcover__{class_group}__within_{buffer_distance}m.gpkg",
    )
    if os.path.exists(output_gpkg) and os.path.exists(output_csv):
        print(f"Skipped: {class_group}__within_{buffer_distance}m")
        return

    print(f"Started: {class_group}__within_{buffer_distance}m")

    # Read
    landcover = geopandas.read_file(
        os.path.join(input_path, "2013_landuse_landcover_split_1k.gpkg")
    )
    landcover = landcover[["OBJECTID", "geometry", "Classify", "index_i", "index_j"]]
    landcover.rename(
        columns={"OBJECTID": "landcover_id", "Classify": "landcover_class"},
        inplace=True,
    )
    landcover["area_ha"] = convert_m2_to_ha(landcover.geometry.area)
    cell_indices = landcover.apply(
        lambda row: (row.index_i, row.index_j), axis=1
    ).unique()

    # Calculate
    classes = CLASS_GROUPS[class_group]

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning)
        group_buffer_summary, group_buffer_surroundings = buffer_and_summarise(
            landcover, cell_indices, class_group, classes, buffer_distance
        )

    # Write
    group_buffer_summary.to_csv(output_csv)
    group_buffer_surroundings.to_file(
        output_gpkg,
        driver="GPKG",
    )
    print(f"Done: {class_group}__within_{buffer_distance}m")


def buffer_and_summarise(
    landcover, cell_indices, class_group, classes, buffer_distance
):
    cell_summaries = []
    surroundings = []
    for i, j in tqdm(cell_indices):
        cell_summary, surrounding = buffer_class_and_intersect(
            select_cell(landcover, i, j),
            select_halo(landcover, i, j),
            classes,
            buffer_distance,
        )
        cell_summaries.append(cell_summary)
        surroundings.append(surrounding)
    summary = (
        pandas.concat(cell_summaries)
        .reset_index()
        .drop(
            columns=[
                "index_i",
                "index_j",
                "cell_area_total_ha",
                "cell_area_surrounding_total_ha",
            ]
        )
        .groupby("landcover_class")
        .sum()
    )

    summary["area_percentage"] = (summary.area_ha / summary.area_ha.sum()) * 100
    summary["buffer_distance_m"] = buffer_distance
    summary["buffer_from_group"] = class_group

    surroundings = pandas.concat(surroundings)
    surroundings["buffer_distance_m"] = buffer_distance
    surroundings["buffer_from_group"] = class_group
    return summary, surroundings


def collection_to_multipolygon(collection):
    return shapely.MultiPolygon(
        geom
        for geom in collection.geoms
        if geom.geom_type in ("Polygon", "MultiPolygon")
    )


def buffer_class_and_intersect(cell, halo, classes, buffer_distance):
    """Buffer areas with `landcover_class` in `classes` to intersect with
    any areas of other `landcover_class` within `buffer_distance`.

    Cell edge length must be <= `buffer_distance`, so that any patch of
    landcover buffered from anywhere in the halo will certainly reach our cell.

    Parameters
    ----------
    cell: geopandas.GeoDataFrame
        contains areas in the grid cell of interest
    halo: : geopandas.GeoDataFrame
        contains areas in the grid cell and all neighbouring cells
    classes: list[str]
        landcover_class to buffer
    buffer_distance: int|float
        buffer distance in metres

    Returns
    -------
    summary: pandas.DataFrame
        with index ["landcover_class", "index_i", "index_j"] and columns
        ["area_ha", "cell_area_total_ha", "cell_area_surrounding_total_ha"]
        containing the total area of each landcover_class in `cell` within
        `buffer_distance` of any landcover with `classes` (with cell indexes for
        reference).
    surrounding_landcover: geopandas.GeoDataFrame
        containing landcover areas in `cell` within `buffer_distance` of any
        landcover with `classes`
    """

    # find any landcover of `classes` in the cell or surrounding halo - these
    # are the areas we will buffer
    landcover_of_class = halo[halo.landcover_class.isin(classes)].copy()

    # find any landcover not of `classes` in the cell - these are our candidate
    # areas for the "buffer strip"
    landcover_not_of_class = cell[~cell.landcover_class.isin(classes)]

    # if there's no landcover of interest, return early with empty DataFrames
    if len(landcover_of_class) == 0:
        return pandas.DataFrame({}), pandas.DataFrame({})

    # find the buffer strip around the landcover of interest
    landcover_geom = landcover_of_class.geometry.unary_union
    buffered_geom = landcover_geom.buffer(distance=buffer_distance)
    surrounding_geom = shapely.difference(buffered_geom, landcover_geom, grid_size=1)

    # We're only interested in the Polygon/MultiPolygon buffer strip - ignore
    # any Point/Linestring features resulting from buffer/difference operation
    if surrounding_geom.geom_type == "GeometryCollection":
        surrounding_geom = shapely.MultiPolygon(
            geom for geom in surrounding_geom.geoms if geom.geom_type == "Polygon"
        )

    # find the other landcover areas within the buffer strip
    surrounding_landcover = landcover_not_of_class.copy()
    surrounding_landcover.geometry = shapely.intersection(
        surrounding_geom, landcover_not_of_class.geometry, grid_size=1
    )

    # Again, we're only interested in Polygon and MultiPolygon results of
    # intersection.
    #
    # For reference see geopandas.tools.overlay with keep_geom_type=True see
    # https://github.com/geopandas/geopandas/blob/d5add48e966bb4e711d2399721df88349a358904/geopandas/tools/overlay.py#L329-L392
    is_collection = surrounding_landcover.geometry.geom_type == "GeometryCollection"
    if is_collection.any():
        collections = surrounding_landcover.geometry[is_collection]
        multipolygons = collections.apply(collection_to_multipolygon)
        surrounding_landcover.loc[is_collection, "geometry"] = multipolygons.values

    surrounding_landcover = surrounding_landcover[
        surrounding_landcover.geometry.geom_type.isin(("Polygon", "MultiPolygon"))
    ]

    # Calculate area in hectares
    surrounding_landcover["area_ha"] = convert_m2_to_ha(
        surrounding_landcover.geometry.area
    )

    # Summarise area statistics
    summary = (
        surrounding_landcover[["landcover_class", "area_ha", "index_i", "index_j"]]
        .groupby(["landcover_class", "index_i", "index_j"])
        .sum()
    )
    summary["cell_area_total_ha"] = cell.area_ha.sum()
    summary["cell_area_surrounding_total_ha"] = summary.area_ha.sum()

    return summary, surrounding_landcover


def convert_m2_to_ha(area):
    return area * 0.0001


def select_cell(df, i, j):
    return df[(df.index_i == i) & (df.index_j == j)]


def select_halo(df, i, j):
    return df[
        ((df.index_i == i - 1) & (df.index_j == j - 1))
        | ((df.index_i == i - 1) & (df.index_j == j))
        | ((df.index_i == i - 1) & (df.index_j == j + 1))
        | ((df.index_i == i) & (df.index_j == j - 1))
        | ((df.index_i == i) & (df.index_j == j))
        | ((df.index_i == i) & (df.index_j == j + 1))
        | ((df.index_i == i + 1) & (df.index_j == j - 1))
        | ((df.index_i == i + 1) & (df.index_j == j))
        | ((df.index_i == i + 1) & (df.index_j == j + 1))
    ]


if __name__ == "__main__":
    try:
        input_path = sys.argv[1]
        output_path = sys.argv[2]
        class_group = sys.argv[3]
        buffer_distance = int(sys.argv[4])
    except IndexError:
        print(
            f"Usage: python {__file__} <input_path> <output_path> <class_group> <buffer_distance>"
        )
        sys.exit(-1)
    main(input_path, output_path, class_group, buffer_distance)
