"""
    Allocating GDP to various asset layers including
"""

rule population_growth_employment:
    """
    Create population projections for Jamaica
    """
    input:
        script = "workflow/1_context/socio_econ_model_scripts/population_growth_employment.py",
        admin_boundaries = f"{DATA}/boundaries/admin_boundaries.gpkg",
        population = f"{DATA}/population/population.gpkg",
        parish_population = f"{RAW}/macroeconomic_data/parish_population_changes.xlsx",
        population_employment = f"{RAW}/macroeconomic_data/population_employment_percent_by_age.csv"
    params:
        epsg = config["adaptation_options"]["epsg_jamaica"],
        base_year = config["economics"]["base_year"]
    output:
        population_projections = f"{DATA}/population/population_projections.gpkg"
    shell:
        """
        python {input.script} \
            --incoming-data-dir {RAW} \
            --data-dir {DATA} \
            --epsg {params.epsg} \
            --base-year {params.base_year} \
            --admin-boundaries {input.admin_boundaries} \
            --population {input.population} \
            --parish-population {input.parish_population} \
            --population-employment {input.population_employment}
        """

rule agriculture_areas_economic_activity:
    """
    Assign agriculture GDP to land use layers in Jamaica
    """
    input:
        script = "workflow/1_context/socio_econ_model_scripts/agriculture_areas_economic_activity.py",
        agri_crop_details = f"{DATA}/agriculture_data/crop_details.csv",
        econ_output = f"{DATA}/macroeconomic_data/detailed_sector_GVA_GDP_current_prices.xlsx",
        land_use = f"{DATA}/land_type_and_use/jamaica_land_use_combined_with_sectors.gpkg",
        intermediate_file = f"{RAW}/buildings/buildings_assigned_economic_sectors_intermediate.gpkg",
        agri_data_prod_dir = f"{RAW}/agriculture_data/spam2010v2r0_global_val_prod_agg.geotiff/JAM",
        agri_data_prod_agg_dir = f"{RAW}/agriculture_data/spam2010v2r0_global_prod.geotiff/JAM",
        agri_data_yield_dir = f"{RAW}/agriculture_data/spam2010v2r0_global_yield.geotiff/JAM",
    params:
        epsg = config["adaptation_options"]["epsg_jamaica"],
        financial_year = config["economics"]["financial_year"]
    output:
        spam_agri_outputs = f"{DATA}/agriculture_data/spam_agriculture_outputs.gpkg",
        prod_keys = f"{DATA}/agriculture_data/production_column_keys.csv",
        tonnage_keys = f"{DATA}/agriculture_data/tonnage_column_keys.csv",
        yield_keys = f"{DATA}/agriculture_data/yield_column_keys.csv",
        agri_gdp = f"{DATA}/agriculture_data/agriculture_gdp.gpkg",
        building_gdp = f"{DATA}/agriculture_data/building_agriculture_gdp.csv"
    shell:
        """
        python {input.script} \
            --data-dir {DATA} \
            --epsg {params.epsg} \
            --financial-year {params.financial_year} \
            --agri-crop-details {input.agri_crop_details} \
            --econ-output {input.econ_output} \
            --land-use {input.land_use} \
            --intermediate-file {input.intermediate_file} \
            --agri-data-prod-dir {input.agri_data_prod_dir} \
            --agri-data-prod-agg-dir {input.agri_data_prod_agg_dir} \
            --agri-data-yield-dir {input.agri_data_yield_dir}
        """

rule mining_areas_economic_activity:
    """
    Assign mining GDP to land use layers in Jamaica
    """
    input:
        script = "workflow/1_context/socio_econ_model_scripts/mining_areas_economic_activity.py",
        mining_gdp = f"{DATA}/mining_data/mining_gdp.gpkg",
        intermediate_file = f"{RAW}/buildings/buildings_assigned_economic_sectors_intermediate.gpkg",
    params:
        epsg = config["adaptation_options"]["epsg_jamaica"],
    output:
        mining_gdp_output = f"{DATA}/mining_data/mining_gdp_with_buildings.gpkg", 
        #originially the above output was the smae mining_gdp.gpk and the name has bee changed
        #An issue arises as to what scripts require this modifie minign gdp as an input as oposed which scripts need mining_gdp before the chanegs 
        building_mining_gdp = f"{DATA}/mining_data/building_mining_gdp.csv"
    shell:
        """
        python {input.script} \
            --epsg {params.epsg} \
            --mining-gdp {input.mining_gdp} \
            --intermediate-file {input.intermediate_file} \
            --mining-gdp-output {output.mining_gdp_output} \
            --building-mining-gdp {output.building_mining_gdp}
        """