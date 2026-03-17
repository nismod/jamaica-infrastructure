"""
Electricity network preprocessing and flow generation.
"""


rule node_source_sink_magnitudes:
    """
    Gather nodal generation and consumption data for electricity network.
    
    Test with:
    snakemake -c1 processed_data/networks/energy/generated_nodal_flows.csv
    """
    input:
        script = "workflow/3b_electricity/generate_nodal_flows.py",
        network = f"{DATA}/networks/energy/electricity_network_v3.2.gpkg",
    output:
        flows = f"{DATA}/networks/energy/generated_nodal_flows.csv",
    shell:
        """
        python {input.script} \\
            --network-file {input.network} \\
            --output-file {output.flows}
        """
