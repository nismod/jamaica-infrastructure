"""Generate hazard-damage curves
"""
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
from matplotlib.ticker import (MaxNLocator,LinearLocator, MultipleLocator)
import matplotlib.pyplot as plt
from matplotlib import cm
from plot_utils import *
from tqdm import tqdm
tqdm.pandas()

mpl.style.use('ggplot')
mpl.rcParams['font.size'] = 10.
mpl.rcParams['font.family'] = 'tahoma'
mpl.rcParams['axes.labelsize'] = 12.
mpl.rcParams['xtick.labelsize'] = 15.
mpl.rcParams['ytick.labelsize'] = 15.

def convert_to_usd(x,loss_column):
    if ("$J" in str(x.damage_cost_unit)) or ("J$" in str(x.damage_cost_unit)) or ("JD" in str(x.damage_cost_unit)):
        return jamaica_currency_conversion()*x[loss_column]
    else:
        return x[loss_column]

def convert_to_jd(x,loss_column):
    if ("$US" in str(x.damage_cost_unit)) or ("US$" in str(x.damage_cost_unit)) or ("USD" in str(x.damage_cost_unit)):
        return (1.0/jamaica_currency_conversion())*x[loss_column]
    else:
        return x[loss_column]


def get_damage_columns(damage_dictionary,damage_dataframe):
    if damage_dictionary["damage_type"] == "direct_damages":
        damage_columns = ["exposure"]
    else:
        damage_columns = []
    for dg_c in damage_dictionary["damage_columns"]:
        damage_columns += [c for c in damage_dataframe.columns.values.tolist() if dg_c in c]
    return damage_columns

def damages_grouped(damages,damage_groupby,damage_sum_columns,sector_name=None,layer_name=None):
    damages = damages.groupby(
                    damage_groupby,dropna=False
                    ).agg(
                        dict(
                            zip(
                                damage_sum_columns,["sum"]*len(damage_sum_columns)
                                )
                            )
                        ).reset_index() 
    damages["sector"] = sector_name
    damages["subsector"] = layer_name
    return damages

def main(config):
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    output_data_path = config['paths']['output']
    figures_data_path = config['paths']['figures']

    hazards = ["fluvial","surface","coastal","cyclone"]
    hazard_labels = ["Fluvial flooding","Pluvial flooding","Coastal flooding","Cyclone winds"]
    hazard_colors = ["#9ecae1","#9e9ac8","#02818a","#78c679"]

    damage_data_path = os.path.join(output_data_path,
                        "direct_damages_summary")

    # Jamaica GDP from 2020, Source: https://statinja.gov.jm/NationalAccounting/Quarterly/NewQuarterlyGDP.aspx
    jamaica_gdp = 2121087*1e6 
    ead_results = pd.read_csv(os.path.join(damage_data_path,
                        "hazard_rcp_epoch_EAD_without_confidence_value_JD.csv"))

    ead_columns = [c for c in ead_results.columns.values.tolist() if "EAD_" in c]
    ear_columns = [c for c in ead_results.columns.values.tolist() if "EAR_" in c]

    total_baseline_eads = ead_results[ead_results["epoch"] == 2010][ead_columns].sum().reset_index()
    total_baseline_eads.columns = ["EAD","values"]
    total_baseline_eads["value_billions"] = 1.0e-9*total_baseline_eads["values"]
    total_baseline_eads["percentage_gdp"] = 100.0*total_baseline_eads["values"]/jamaica_gdp
    # print (total_baseline_eads)

    total_baseline_ears = ead_results[ead_results["epoch"] == 2010][ear_columns].sum().reset_index()
    total_baseline_ears.columns = ["EAR","values"]
    total_baseline_ears["value_billions"] = 1.0e-9*total_baseline_ears["values"]
    total_baseline_ears["percentage_gdp"] = 100.0*total_baseline_ears["values"]/jamaica_gdp
    # print (total_baseline_ears)
    hazard_df = []
    for hazard in hazards:
        hazard_eads = ead_results[(ead_results["epoch"] == 2010) & (ead_results["hazard"] == hazard)][ead_columns].sum().reset_index()
        hazard_eads.columns = ["EAD",hazard]
        if len(hazard_df) > 0:
            hazard_df = pd.merge(hazard_df,hazard_eads,how="left",on=["EAD"])
        else:
            hazard_df = hazard_eads.copy()

    # print (hazard_df.head(3))

    # hazard_ears = ead_results[ead_results["epoch"] == 2010].groupby("hazard")[ear_columns].sum().reset_index()
    # print (hazard_ears)

    fig, ax = plt.subplots(1,1,figsize=(8,8),dpi=500)
    boxp = ax.boxplot(1e-9*hazard_df[hazards].head(3),showfliers=False,patch_artist=True,labels=hazards)
    for patch, color in zip(boxp['boxes'], hazard_colors):
        patch.set_facecolor(color)
    plt.setp(boxp['medians'], color="#000000")
    
    ax.set_ylabel('Expected Annual Damages (J$ Billion)',fontweight='bold',fontsize=15)
    ax.set_xticklabels(hazard_labels,fontsize=15, rotation=0,fontweight="bold")
    plt.tight_layout()
    save_fig(os.path.join(figures_data_path,
                'jamaica_hazard_EAD_totals.png'))
    plt.close() 

    # print (total_baseline_ears)
    sectors = list(set(ead_results["sector"].values.tolist()))
    hazard_df = []
    # sector_labels = []
    for sector in sectors:
        hazard_eads = ead_results[(ead_results["epoch"] == 2010) & (ead_results["sector"] == sector)][ead_columns].sum().reset_index()
        hazard_eads.columns = ["EAD",sector]
        if len(hazard_df) > 0:
            hazard_df = pd.merge(hazard_df,hazard_eads,how="left",on=["EAD"])
        else:
            hazard_df = hazard_eads.copy()

    # print (hazard_df.head(3))

    # hazard_ears = ead_results[ead_results["epoch"] == 2010].groupby("hazard")[ear_columns].sum().reset_index()
    # print (hazard_ears)

    fig, ax = plt.subplots(1,1,figsize=(8,8),dpi=500)
    boxp = ax.boxplot(1e-9*hazard_df[sectors].head(3),showfliers=False,patch_artist=True,labels=sectors,vert=False)
    # for patch, color in zip(boxp['boxes'], hazard_colors):
    #     patch.set_facecolor(color)
    plt.setp(boxp['medians'], color="#000000")
    
    ax.set_xlabel('Expected Annual Damages (J$ Billion)',fontweight='bold',fontsize=15)
    ax.set_yticklabels(sectors,fontsize=15, rotation=0,fontweight="bold")
    plt.tight_layout()
    save_fig(os.path.join(figures_data_path,
                'jamaica_hazard_sector_EAD_totals.png'))
    plt.close() 

    # print (total_baseline_ears)
    for i, (hazard,hazard_color,hazard_label)  in enumerate(list(zip(hazards,hazard_colors,hazard_labels))):
        sectors = list(set(ead_results["sector"].values.tolist()))
        hazard_df = []
        # sector_labels = []
        for sector in sectors:
            hazard_eads = ead_results[
                                    (
                                        ead_results["epoch"] == 2010
                                    ) & (
                                        ead_results["hazard"] == hazard
                                        ) & (
                                            ead_results["sector"] == sector
                                            )][ead_columns].sum().reset_index()
            hazard_eads.columns = ["EAD",sector]
            if len(hazard_df) > 0:
                hazard_df = pd.merge(hazard_df,hazard_eads,how="left",on=["EAD"])
            else:
                hazard_df = hazard_eads.copy()

        # print (hazard_df.head(3))

        # hazard_ears = ead_results[ead_results["epoch"] == 2010].groupby("hazard")[ear_columns].sum().reset_index()
        # print (hazard_ears)

        fig, ax = plt.subplots(1,1,figsize=(8,8),dpi=500)
        boxp = ax.boxplot(1e-9*hazard_df[sectors].head(3),showfliers=False,patch_artist=True,labels=sectors,vert=False)
        for patch, color in zip(boxp['boxes'], [hazard_color]*len(boxp['boxes'])):
            patch.set_facecolor(color)
        plt.setp(boxp['medians'], color="#000000")
        
        ax.set_xlabel('Expected Annual Damages (J$ Billion)',fontweight='bold',fontsize=15)
        ax.set_yticklabels(sectors,fontsize=15, rotation=0,fontweight="bold")
        ax.set_title(f"{hazard_label} Physical Risks",fontweight="bold",fontsize=15)
        plt.tight_layout()
        save_fig(os.path.join(figures_data_path,
                    f'jamaica_hazard_sector_EAD_totals_{hazard}.png'))
        plt.close()

    damage_results = pd.read_csv(os.path.join(damage_data_path,
                        "hazard_rcp_epoch_direct_damages_without_confidence_value_JD.csv"))

    damage_columns = [c for c in damage_results.columns.values.tolist() if "direct_damage_cost" in c]
    reopen_columns = [c for c in damage_results.columns.values.tolist() if "direct_reopen_cost" in c]

    total_baseline_damages = damage_results[(damage_results["epoch"] == 2010) & (damage_results["rp"] == 100)][damage_columns].sum().reset_index()
    total_baseline_damages.columns = ["damage_cost","values"]
    total_baseline_damages["value_billions"] = 1.0e-9*total_baseline_damages["values"]
    total_baseline_damages["percentage_gdp"] = 100.0*total_baseline_damages["values"]/jamaica_gdp
    # print (total_baseline_damages)

    hazard_df = []
    for hazard in hazards:
        hazard_damages = damage_results[
                            (
                                damage_results["epoch"] == 2010
                                ) & (
                                    damage_results["hazard"] == hazard
                                    ) & (
                                        damage_results["rp"] == 100
                                        )][damage_columns].sum().reset_index()
        hazard_damages.columns = ["damage_cost",hazard]
        if len(hazard_df) > 0:
            hazard_df = pd.merge(hazard_df,hazard_damages,how="left",on=["damage_cost"])
        else:
            hazard_df = hazard_damages.copy()

    # print (hazard_df.head(3))

    # hazard_ears = ead_results[ead_results["epoch"] == 2010].groupby("hazard")[ear_columns].sum().reset_index()
    # print (hazard_ears)
    rp = 100
    fig, ax = plt.subplots(1,1,figsize=(8,8),dpi=500)
    boxp = ax.boxplot(1e-9*hazard_df[hazards].head(3),showfliers=False,patch_artist=True,labels=hazards)
    for patch, color in zip(boxp['boxes'], hazard_colors):
        patch.set_facecolor(color)
    plt.setp(boxp['medians'], color="#000000")
    
    ax.set_ylabel(f'Total Damages for 1 in {rp}-year event (J$ Billion)',fontweight='bold',fontsize=15)
    ax.set_xticklabels(hazard_labels,fontsize=15, rotation=0,fontweight="bold")
    plt.tight_layout()
    save_fig(os.path.join(figures_data_path,
                f'jamaica_hazard_1in{rp}_event_totals.png'))
    plt.close()

    hazard_colors = ["#08519c","#54278f","#016c59","#006d2c"]
    for i, (hazard, hazard_label,hazard_color) in enumerate(list(zip(hazards,hazard_labels,hazard_colors))):
        hazard_df = damage_results[
                            (
                                damage_results["epoch"] == 2010
                                ) & (
                                    damage_results["hazard"] == hazard
                                    )].groupby(['hazard','rp'])[damage_columns].sum().reset_index()
        # print (hazards_df)

        fig, ax = plt.subplots(1,1,figsize=(8,8),dpi=500)
        ax.plot(hazard_df['rp'],
                    1e-9*hazard_df["direct_damage_cost_mean"],'-',
                    marker='o',
                    color=hazard_color,linewidth=2.0,
                    label=f'{hazard_label} damages')    

        ax.fill_between(hazard_df['rp'],
                        1e-9*hazard_df["direct_damage_cost_min"],
                        1e-9*hazard_df["direct_damage_cost_max"],
                        alpha=0.5,facecolor=hazard_color,
                        label='Uncertainty range')

        ax.legend(loc='upper left',fontsize=14)
        ax.set_xlabel("Return period (years)",fontweight='bold',fontsize=12)
        ax.set_ylabel('Direct damages (J$ billion)',fontweight='bold',fontsize=12)

        plt.tight_layout()
        save_fig(os.path.join(figures_data_path,
                f'jamaica_return_period_damages_{hazard}.png'))
        plt.close()


    for i, (hazard, hazard_label,hazard_color) in enumerate(list(zip(hazards,hazard_labels,hazard_colors))):
        hazard_df = damage_results[
                            (
                                damage_results["epoch"] == 2010
                                ) & (
                                    damage_results["hazard"] == hazard
                                    )].groupby(['hazard','rp'])[damage_columns].sum().reset_index()
        # print (hazards_df)
        if hazard == 'cyclone':
            fig, ax = plt.subplots(1,1,figsize=(8,8),dpi=500)
            # ax.plot(1.0/hazard_df['rp'],
            #             1e-9*hazard_df["direct_damage_cost_mean"],'-',
            #             marker='o',
            #             color="#cb181d",linewidth=2.0,
            #             label=f'{hazard_label} damages')
            ax.plot(1.0/hazard_df['rp'],
                        1e-9*hazard_df["direct_damage_cost_mean"],'-',
                        marker='o',
                        color="#000000",linewidth=2.0,
                        label="Damage or Loss vs Exceedance probability curve")    

            ax.fill_between(1.0/hazard_df['rp'],
                            np.array([0]*len(hazard_df.index)),
                            1e-9*hazard_df["direct_damage_cost_mean"],
                            alpha=0.5,facecolor="#737373",
                            label="Area under curve")

            ax.legend(loc='upper right',fontsize=14)
            ax.set_xlabel("Exceedance probability",fontweight='bold',fontsize=16)
            # ax.set_ylabel('Direct damages (J$ billion)',fontweight='bold',fontsize=12)
            ax.set_ylabel('Damages or Losses',fontweight='bold',fontsize=16)

            plt.tight_layout()
            # save_fig(os.path.join(figures_data_path,
            #         f'jamaica_probability_damages_{hazard}.png'))
            save_fig(os.path.join(figures_data_path,
                    "example_loss_probability_curve.png"))
            plt.close()


    """Climate change plots
    """
    currency_inflation = 0.05
    rcps = [2.6,4.5,8.5]
    baseline_year = 2010
    rcp_colors = ["#006d2c","#08519c","#b30000"]
    ead_results.loc[ead_results["epoch"] == 2010,"rcp"] = "baseline"
    for i, (hazard,hazard_color,hazard_label)  in enumerate(list(zip(hazards,hazard_colors,hazard_labels))):
        if hazard in ["cyclone","coastal"]:
            fig, ax = plt.subplots(1,1,figsize=(8,8),dpi=500)
            hazard_df = ead_results[ead_results["hazard"] == hazard]
            hazard_df = hazard_df.groupby(["rcp","epoch"])[ead_columns].sum().reset_index()
            hazard_baseline = hazard_df[hazard_df["rcp"] == "baseline"]
            for j, (rcp,rcp_color) in enumerate(list(zip(rcps,rcp_colors))):
                hazard_rcp = hazard_df[hazard_df["rcp"] == str(rcp)]
                if len(hazard_rcp.index) > 0:
                    hazard_rcp["EAD_undefended_mean"] = hazard_rcp["EAD_undefended_mean"]*(1+currency_inflation)**(hazard_rcp["epoch"] - baseline_year)
                    hazard_rcp = pd.concat([hazard_baseline,hazard_rcp],axis=0,ignore_index=True)
                    ax.plot(hazard_rcp['epoch'],
                            1e-9*hazard_rcp["EAD_undefended_mean"],'-',
                            marker='o',
                            color=rcp_color,linewidth=2.0,
                            label=f'RCP {rcp}')    

            ax.legend(loc='upper left',fontsize=14)
            ax.set_xlabel("Time epoch (years)",fontweight='bold',fontsize=12)
            ax.set_ylabel('Average Expected Annual Damages (J$ billion)',fontweight='bold',fontsize=12)

            plt.tight_layout()
            save_fig(os.path.join(figures_data_path,
                    f'jamaica_average_ead_{hazard}_climate_scenarios_inflation.png'))
            plt.close()



if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
