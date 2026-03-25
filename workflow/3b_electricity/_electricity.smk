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
        network = f"{DATA}/networks/energy/electricity_network_{{version}}.gpkg",
    output:
        # TODO: Criticality and hotspot rules use this as input, but do not currently discriminate version
        flows = f"{DATA}/networks/energy/generated_nodal_flows_{{version}}.csv",
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
    snakemake -c1 results/electricity_failures/diagnostics/edge_flows_v3.2.geoparquet
    """
    input:
        script = "workflow/3b_electricity/nominal_flow.py",
        network = f"{DATA}/networks/energy/electricity_network_{{version}}.gpkg",
        flows = f"{DATA}/networks/energy/generated_nodal_flows_{{version}}.csv",
    output:
        edges = f"{OUTPUT}/electricity_failures/diagnostics/edges_{{version}}.geoparquet",
        nodes = f"{OUTPUT}/electricity_failures/diagnostics/nodes_{{version}}.geoparquet",
    shell:
        """
        python {input.script} \\
            --network-path {input.network} \\
            --flows-path {input.flows} \\
            --output-edges-path {output.edges} \\
            --output-nodes-path {output.nodes}
        """
