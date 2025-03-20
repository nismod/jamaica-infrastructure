"""
- Calculating expected annual damages (EAD) from direct damages 
- Sector specific losses
- Calculating expected annual economic losses (EAEL)
"""


rule compute_EAD_and_EAEL:
    """
    Calculate expected annual damages and expected annual economic losses for all networks and hazards.
    Test with:
    snakemake -c1 results/direct_damages/EAD_EAEL.flag
    """
    input:
        damage_results = f"{OUTPUT}/direct_damages/ead_eael_results.txt",
        script_eal_eael = "scripts/analysis/ead_eael_calculations.py",
        # we need per edge single asset failure costs to calculate EAEL
        # there are precomputed files here:
        # linux-filestore://soge-home/projects/mistral/jamaica-ccri/results/economic_losses/single_failure_scenarios
        # must recreate for any network with:
        # -> scripts/analysis/*_single_point_failure_results_combine.py
        # depends on:
        # -> scripts/transport_model/transport_failure_scenario_setup.py
        # single_point_failure = ?
    output:
        eal_eael_flag = f"{OUTPUT}/direct_damages/EAD_EAEL.flag"
    shell:
        """
        parallel --lb -j {workflow.cores} --colsep "," -a {input.damage_results} python {input.script_eal_eael} {{}}
        """


rule summarise_damages_and_losses:
    """
    Summarise direct damages and losses.
    Test with:
    snakemake -c1 results/direct_damages_summary/summary.flag
    """
    input:
        script_loss_summary = "scripts/analysis/damage_loss_summarised.py",
        eal_eael_flag = f"{OUTPUT}/direct_damages/EAD_EAEL.flag",
        damage_results = f"{OUTPUT}/direct_damages",
        networks = config["paths"]["network_layers"],
        sensitivity_parameters = f"{DATA}/sensitivity_parameters.csv"
    output:
        summary_flag = f"{OUTPUT}/direct_damages_summary/summary.flag"
    shell:
        f"""
        python {{input.script_loss_summary}} {{input.damage_results}} {OUTPUT}/direct_damages_summary \
            {{input.networks}} {{input.sensitivity_parameters}}
        """


rule calculate_NPV:
    """
    Calculate the Net Present Value and timeseries of risks, given discounting
    and projected GDP growth.
    Test with:
    snakemake -c1 results/direct_damages_summary/losses_npv.flag
    """
    input:
        script_loss_npv = "scripts/analysis/damage_loss_timeseries_and_npv.py",
    output:
        loss_npv_flag = f"{OUTPUT}/direct_damages_summary/losses_npv.flag"
    shell:
        f"""
        python {{input.script_loss_npv}} {OUTPUT}/direct_damages_summary {OUTPUT}/loss_damage_timeseries \
            {OUTPUT}/loss_damage_npvs {config["paths"]["network_layers"]}
        """
