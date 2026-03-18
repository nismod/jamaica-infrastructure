"""Prepare building footprints with attributes for allocation to economic sectors.
"""

rule download_osm_jamaica:
    """Download country extract
    """
    output:
        pbf = f"{RAW}/osm/jamaica-{config["context"]["osm_extract_date"]}.osm.pbf",
    shell:
        """
        extract_dir=$(dirname {output.pbf})
        extract_date={config[context][osm_extract_date]}

        pushd $extract_dir

            # download extract
            wget http://download.geofabrik.de/central-america/jamaica-${{extract_date}}.osm.pbf
            wget http://download.geofabrik.de/central-america/jamaica-${{extract_date}}.osm.pbf.md5

            # check extract
            md5sum --check jamaica-${{extract_date}}.osm.pbf.md5

        popd
        """


rule filter_osm_buildings:
    """Filter OpenStreetMap extract for all buildings in Jamaica
    """
    input:
        pbf = f"{RAW}/osm/jamaica-{config["context"]["osm_extract_date"]}.osm.pbf",
    output:
        pbf = f"{DATA}/buildings/jamaica-buildings.osm.pbf",
        gpkg = f"{DATA}/buildings/jamaica-buildings.gpkg",
    shell:
        """
        osmium tags-filter \
            {input.pbf} \
            wnr/building \
            --overwrite \
            -o {output.pbf}

        # Extract features from .osm.pbf to .gpkg
        OSM_CONFIG_FILE=workflow/1_context/osmconf.ini ogr2ogr -f GPKG \
            {output.gpkg} \
            {output.pbf} \
            multipolygons
        """


rule parse_nic_2005_building_sectors:
    input:
        script = "workflow/1_context/parse_nic_2005_building_sectors.py",
        buildings = f"{DATA}/buildings/buildings_assigned_economic_activity.archive-202603.gpkg",
    output:
        buildings = f"{DATA}/buildings/buildings_nic2005.geoparquet",
    shell:
        """
        python {input.script} \
            --buildings {input.buildings} \
            --output {output.buildings}
        """


rule tag_buildings_with_nic_2016_sector_code:
    """
    Given a set of two digit Jamaica Industrial Classification (JIC) 2005
    codes, identify the relevant JIC 2016 two digit code and sector letter.

    Test with:
    snakemake -c1 processed_data/buildings/buildings_nic2016.geoparquet
    """
    input:
        buildings = rules.parse_nic_2005_building_sectors.output.buildings,
        jic_map = f"{DATA}/macroeconomic_data/JIC_2005-2016_conversion_table.csv",
    output:
        buildings = f"{DATA}/buildings/buildings_nic2016.geoparquet"
    run:
        import geopandas as gpd
        import pandas as pd

        jic_map = pd.read_csv(input.jic_map, comment="#", dtype=str).loc[:, ["jic2016_l1", "jic2016_l2", "jic2005_l2"]]
        lookup = dict(zip(jic_map.jic2005_l2, jic_map.jic2016_l1))

        build = gpd.read_parquet(input.buildings)
        build = build.reset_index(drop=True)

        df = build.loc[:, ["osm_id", "jic2005_two_digit"]].copy()
        df["jic2005_two_digit_list"] = df.jic2005_two_digit.apply(lambda s: list(s.split(",")))
        df = df.explode("jic2005_two_digit_list")
        df = df.rename(columns={"jic2005_two_digit_list": "jic2005_l2"})
        df["jic2016_sector"] = df["jic2005_l2"].map(lookup)
        df = df.dropna(subset="jic2016_sector")
        df = df.loc[:, ["osm_id", "jic2016_sector"]].groupby("osm_id").agg(lambda x: ",".join(x))

        join = build.merge(df.reset_index(), on="osm_id", how="outer")
        join.to_parquet(output.buildings)

