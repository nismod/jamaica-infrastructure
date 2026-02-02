import pandas as pd


def economic_losses_from_network_damage(
    results: pd.DataFrame,
    id_col: str,
    no_access_col: str,
    all_flows: pd.DataFrame,
    trade_sectors: list[str],
    trade_effect: float,
    hourly_wage: float,
) -> pd.DataFrame:

    index_cols: list[str] = [id_col, no_access_col]

    results = pd.merge(results, all_flows, how="left", on=["origin_id", "destination_id"]).fillna(0)
    results["total_trade"] = results[[f"{t}_trade" for t in trade_sectors]].sum(axis=1)

    results["time_loss"] = (1 - results[no_access_col]) * (results["new_cost"] - results["gcost"])
    results["labour_rerouting_loss"] = hourly_wage * results["time_loss"] * results["working_trips"]
    results["trade_rerouting_loss"] = trade_effect * results["time_loss"] * results["total_trade"]
    results["labour_gdp_loss"] = results[no_access_col] * results["GDP_to_trips"]
    results["trade_loss"] = results[no_access_col] * results["total_trade"]

    losses = (
        results.loc[
            :,
            index_cols
            + [
                "time_loss",
                "labour_rerouting_loss",
                "trade_rerouting_loss",
                "labour_gdp_loss",
                "trade_loss",
            ],
        ]
        .groupby(index_cols)
        .sum()
        .reset_index()
    )

    rerouting_times_min = results.loc[:, index_cols + ["time_loss"]].groupby(index_cols).min().reset_index()
    rerouting_times_min = rerouting_times_min.rename(columns={"time_loss": "min_trip_time_loss"})
    rerouting_times_max = results.loc[:, index_cols + ["time_loss"]].groupby(index_cols).max().reset_index()
    rerouting_times_max = rerouting_times_max.rename(columns={"time_loss": "max_trip_time_loss"})
    rerouting_times_mean = results.loc[:, index_cols + ["time_loss"]].groupby(index_cols).mean().reset_index()
    rerouting_times_mean = rerouting_times_mean.rename(columns={"time_loss": "mean_trip_time_loss"})

    losses = pd.merge(losses, rerouting_times_min.drop(columns=[no_access_col]), how="left", on=[id_col])
    losses = pd.merge(losses, rerouting_times_max.drop(columns=[no_access_col]), how="left", on=[id_col])
    losses = pd.merge(losses, rerouting_times_mean.drop(columns=[no_access_col]), how="left", on=[id_col])

    return losses
