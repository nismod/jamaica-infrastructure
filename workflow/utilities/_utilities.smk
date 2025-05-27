"""
Helper functions used across rule files.
"""

def get_n_ensemble(wildcards) -> int:
    """
    Wait for ensemble parameters to be generated (checkpoint), then read the
    length of the list.
    """
    filepath = checkpoints.sensitivity_parameters.get(**wildcards).output.sensitivity_parameters
    return len(pd.read_csv(filepath))

def sensitivity_id_from_slug(wildcards):
    """
    'parameter_set_0' -> '0'
    """
    return wildcards.parameter_set.replace("parameter_set_", "")

def get_asset_metadata(wildcards) -> pd.Series:
    """
    With `asset_gpkg`, e.g. 'port_polygon'  and `asset_layer`, e.g. 'areas'
    wildcards, get the metadata for the asset class from file.
    """
    df = pd.read_csv(config["paths"]["network_layers"])
    row = df[(df['asset_gpkg'] == wildcards.gpkg) & (df['asset_layer'] == wildcards.layer)]
    if len(row) > 1:
        raise ValueError(f"Multiple assets found for gpkg={wildcards.gpkg} and layer={wildcards.layer}")
    elif len(row) == 0:
        raise ValueError(f"No asset found for gpkg={wildcards.gpkg} and layer={wildcards.layer}")
    else:
        return row.squeeze()

def get_flood_threshold_from_path(wildcards) -> float:
    output_path = wildcards.output_path
    match = re.search(r"results/flood_threshold_(\d+)p(\d+)", output_path)
    if match:
        w, d = match.groups()
        return float(f"{w}.{d}")
    else:
        return 0.5
    