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


rule electricity_nominal_flow:
    """
    Runs JEM model with no asset failures and saves flow results to geoparquet
    for visualization and analysis.

    N.B. You will need a valid Gurobi license file. The environmental variable
    GRB_LICENSE_FILE can be used to specify a gurobi.lic file path.
    
    Test with:
    snakemake -c1 results/electricity_failures/diagnostics/nominal_flows.geoparquet
    """
    input:
        script = "workflow/3b_electricity/nominal_flow.py",
        network = f"{DATA}/networks/energy/electricity_network_v3.2.gpkg",
        flows = f"{DATA}/networks/energy/generated_nodal_flows.csv",
    output:
        flows_data = f"{OUTPUT}/electricity_failures/diagnostics/nominal_flows.geoparquet",
        nodes_data = f"{OUTPUT}/electricity_failures/diagnostics/nominal_nodes.geoparquet",
    shell:
        """
        python {input.script} \\
            --network-path {input.network} \\
            --flows-path {input.flows} \\
            --output-path {output.flows_data}
        """
