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