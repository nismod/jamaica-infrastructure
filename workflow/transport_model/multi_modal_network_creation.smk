rule create_multi_modal_network:
    """
    Create a multi-modal transport network.
    Test with:
    snakemake -c1 processed_data/networks/transport/multi_modal_network.gpkg
    """
    input:
        airports = "{PROCESSED_DATA}/networks/transport/airport_polygon.gpkg",
        ports = "{PROCESSED_DATA}/networks/transport/port_polygon.gpkg",
        rail = "{PROCESSED_DATA}/networks/transport/rail.gpkg",
        road = "{PROCESSED_DATA}/networks/transport/roads.gpkg",
    output:
        multi_modal_network = "{PROCESSED_DATA}/networks/transport/multi_modal_network.gpkg",
    script:
        "./multi_modal_network_creation.py"
