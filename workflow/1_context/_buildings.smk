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

        jic_map = pd.read_csv(input.jic_map, comment="#", dtype=str).loc[:, ["jic2016_l1", "jic2005_l2"]]
        # jic005_l2 can map to multiple jic2016_l1
        jic_map = jic_map.groupby("jic2005_l2").agg(lambda cs: ",".join(sorted(set(cs)))).reset_index()

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

rule region_attractiveness:
    """Given regional population, estimate an attractiveness term based on
    population within a distance of each region
    """
    input:
        script = "workflow/1_context/region_attractiveness.py",
        population = f"{DATA}/population/population_projections.gpkg",
    output:
        regions = f"{DATA}/population/region_attractiveness.geoparquet"
    shell:
        """
        python {input.script} \
            --population_path {input.population} \
            --output_region_attractiveness_path {output.regions}
        """


rule agriculture_crops:
    """Process agricultural crop production to areas, from directory of SPAM rasters
    """
    input:
        script = "workflow/1_context/agricultural_crops.py",
        spam_path = f"{RAW}/agriculture_data",
        spam_agriculture_outputs = f"{DATA}/agriculture_data/spam_agriculture_outputs.gpkg",
    shell:
        """
        python {input.script} \
            --spam-path {input.spam_path} \
            --spam-agriculture-outputs-path {input.spam_agriculture_outputs}
        """


rule agriculture_gva:
    """Given agricultural land-use, modelled crop yields, and national GVA for
    agriculture and forestry subsectors, estimate agricultural area production
    values.
    """
    input:
        script = "workflow/1_context/agricultural_production.py",
        spam_agriculture_outputs = f"{DATA}/agriculture_data/spam_agriculture_outputs.gpkg",
        crop_details = f"{DATA}/agriculture_data/crop_details.NIP_2023_codes.csv",
        land_use = f"{DATA}/land_type_and_use/jamaica_land_use_combined_with_sectors.gpkg",
        fishing = f"{DATA}/land_type_and_use/aqua_farms.gpkg",
        economic_output = f"{DATA}/macroeconomic_data/NIP_2023.csv",
    output:
        agriculture_gdp = f"{DATA}/agriculture_data/agriculture_gdp.gpkg"
    shell:
        """
        python {input.script} \
            --spam-agriculture-outputs-path {input.spam_agriculture_outputs} \
            --crop-details-path {input.crop_details} \
            --land-use-path {input.land_use} \
            --fishing-locations-path {input.fishing} \
            --economic-output-path {input.economic_output} \
            --output-areas {output.agriculture_gdp}
        """


rule mining_gva:
    """Given mining land-use and national GVA for quarrying and bauxite extraction,
    estimate mining area production values, weighted by:
    - trade volume of nearest port over land transport (road/rail) network
    - footprint of mining area
    """
    input:
        script = "workflow/1_context/mine_production.py",
        land_use = f"{DATA}/land_type_and_use/jamaica_land_use_combined_with_sectors.gpkg",
        ports = f"{DATA}/networks/transport/port_polygon.gpkg",
        network = f"{DATA}/networks/transport/multi_modal_network.gpkg",
        economic_output = f"{DATA}/macroeconomic_data/NIP_2023.csv",
    output:
        mining_gdp = f"{DATA}/mining_data/mining_gdp.gpkg"
    shell:
        """
        python {input.script} \
            --land_use_path {input.land_use} \
            --ports_path {input.ports} \
            --network_path {input.network} \
            --economic_output_path {input.economic_output} \
            --output_path {output.mining_gdp}
        """

rule allocate_gva:
    """Given buildings tagged by sector, and regional attractivess, estimate
    building daily GDP (JMD/day) - N.B. this is properly GVA, excluding taxes
    and subsidies.

    TODO: handle mining and agriculture reasonably
    """
    input:
        script = "workflow/1_context/spatial_economic_allocation.py",
        buildings = rules.tag_buildings_with_nic_2016_sector_code.output.buildings,
        regions = rules.region_attractiveness.output.regions,
        national_industrial_product = f"{DATA}/macroeconomic_data/NIP_2023.csv",
        agriculture = f"{DATA}/agriculture_data/agriculture_gdp.gpkg",
        mining = f"{DATA}/mining_data/mining_gdp.gpkg",
    output:
        buildings = f"{DATA}/buildings/buildings_assigned_economic_activity.geoparquet"
    shell:
        """
        python {input.script} \
            --buildings_path {DATA}/buildings//buildings_nic2016.geoparquet \
            --region_attractiveness_path {input.regions} \
            --economic_output_path {input.national_industrial_product} \
            --agriculture_areas_path {input.agriculture} \
            --mining_areas_path {input.mining} \
            --output {output.buildings}
        """

rule regional_gva:
    input:
        population = f"{DATA}/population/population_projections.gpkg",
        buildings = rules.allocate_gva.output.buildings
    output:
        regions = f"{DATA}/buildings/admin_level_assigned_economic_activity.geoparquet"
    run:
        import geopandas as gpd
        import pandas as pd
        from jamaica_infrastructure.geo import LOCAL_PROJ_CRS_EPSG

        population_areas = gpd.read_file(input.population, layer="mean")
        buildings = gpd.read_parquet(input.buildings)
        sector_codes = sorted(
            list(
                pd.Series(buildings.jic2016_sector.dropna().unique())
                .apply(lambda s: s.split(","))
                .explode()
                .unique()
            )
        )
        gdp_columns = [f"{scode}_GDP" for scode in sector_codes] + ["total_GDP"]

        admin_gdp = (
            buildings.groupby(["ED_ID", "ED", "PARISH", "CONST_NAME"])[gdp_columns]
            .sum()
            .reset_index()
        )
        admin_gdp["GDP_unit"] = "JD/day"
        admin_gdp = gpd.GeoDataFrame(
            pd.merge(
                admin_gdp,
                population_areas[["ED_ID", "ED", "geometry"]],
                how="left",
                on=["ED_ID", "ED"],
            ),
            geometry="geometry",
            crs=f"EPSG:{LOCAL_PROJ_CRS_EPSG}",
        )

        admin_gdp.to_parquet(output.regions)
