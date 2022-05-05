"""Generate hazard-damage curves
"""
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
# from matplotlib.ticker import (MaxNLocator,LineaelLocator, MultipleLocator)
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

    hazard_plot_df = pd.DataFrame(list(zip(hazards,hazard_labels,hazard_colors)),
                            columns=["hazard","hazard_labels","hazard_colors"])
    
    damage_data_path = os.path.join(output_data_path,
                        "damage_loss_sums")

    # Jamaica GDP from 2019, Source: https://statinja.gov.jm/NationalAccounting/Quarterly/NewQuarterlyGDP.aspx
    jamaica_gdp = 2110433*1e6 
    days = 1
    multiply_factor = 1.0e-9

    ead_results = pd.read_csv(os.path.join(damage_data_path,
                        "hazard_rcp_epoch_EAD_EAEL_totals.csv"))
    # ead_results = ead_results[(ead_results["sector"] != "Potable water") | (ead_results["subsector"] != "edges")]
    print (ead_results)
    
    ead_columns = [c for c in ead_results.columns.values.tolist() if "EAD_" in c]
    eael_columns = [c for c in ead_results.columns.values.tolist() if "EAEL_" in c]
    ead_results[eael_columns] = days*ead_results[eael_columns]

    total_baseline_eads = ead_results[ead_results["epoch"] == 2010][ead_columns].sum().reset_index()
    total_baseline_eads.columns = ["EAD","values"]
    total_baseline_eads["value_billions"] = multiply_factor*total_baseline_eads["values"]
    total_baseline_eads["percentage_gdp"] = 100.0*total_baseline_eads["values"]/jamaica_gdp
    print (total_baseline_eads)

    total_baseline_eaels = ead_results[ead_results["epoch"] == 2010][eael_columns].sum().reset_index()
    total_baseline_eaels.columns = ["eael","values"]
    total_baseline_eaels["value_billions"] = multiply_factor*total_baseline_eaels["values"]
    total_baseline_eaels["percentage_gdp"] = 100.0*total_baseline_eaels["values"]/jamaica_gdp
    print (total_baseline_eaels)
    
    hazard_df = []
    for hazard in hazards:
        hazard_eads = ead_results[(ead_results["epoch"] == 2010) & (ead_results["hazard"] == hazard)][ead_columns].sum().reset_index()
        hazard_eads.columns = ["EAD",hazard]
        if len(hazard_df) > 0:
            hazard_df = pd.merge(hazard_df,hazard_eads,how="left",on=["EAD"])
        else:
            hazard_df = hazard_eads.copy()

    # print (hazard_df.head(3))

    df = ead_results[ead_results["epoch"] == 2010].groupby("hazard")[ead_columns + eael_columns].sum().reset_index()
    # print (df)

    df = pd.merge(hazard_plot_df,df,how="left",on=["hazard"])
    print (df)

    d_color = "#de2d26"
    i_color = "#2c7fb8"
    l_color = "#000000"
    legend_handles = []
    legend_handles.append(mpatches.Patch(color=d_color, label="Direct Damages (Mean Value)"))
    legend_handles.append(mpatches.Patch(color=i_color, label="Indirect Economic Losses (Mean Value)"))
    legend_handles.append(Line2D([0], [0], 
                            color=l_color, label="Min-Max Range"))

    width = 0.35
    fig, ax = plt.subplots(1,1,figsize=(8,8),dpi=500)
    ax.bar(df["hazard"],multiply_factor*df["EAD_mean"],width=width,color=d_color,label="Direct damages")
    ax.bar(df["hazard"],multiply_factor*df["EAEL_mean"],width=width,color=i_color,bottom=multiply_factor*df["EAD_mean"],
                yerr=(multiply_factor*(df["EAD_mean"] - df["EAD_amin"]),multiply_factor*(df["EAD_amax"] - df["EAD_mean"])),
                capsize=5,label="Indirect economic losses")

    # ax.legend(prop={'size':12,'weight':'bold'})
    ax.legend(handles=legend_handles,prop={'size':12,'weight':'bold'})
    ax.set_ylabel('Expected Annual Damages and Losses (J$ Billion)',fontweight='bold',fontsize=15)
    ax.set_xticklabels(df["hazard_labels"].values.tolist(),fontsize=15, rotation=0,fontweight="bold")
    plt.tight_layout()
    save_fig(os.path.join(figures_data_path,
                'jamaica_hazard_EAD_EAEL_totals.png'))
    plt.close()


    df = ead_results[ead_results["epoch"] == 2010].groupby(["hazard","sector"])[ead_columns + eael_columns].sum().reset_index()
    print (df)
    all_sectors = list(set(df["sector"].values.tolist()))
    max_value = multiply_factor*(df["EAD_amax"].max() + df["EAEL_amax"].max()) + 5

    for idx,(hz,hz_l) in enumerate(list(zip(hazards,hazard_labels))):
        df_hz = df[df["hazard"] == hz]
        add_sectors = []
        for s in all_sectors:
            if s not in df_hz["sector"].values.tolist():
                val_dict = dict(zip(ead_columns + eael_columns,[1]*len(ead_columns + eael_columns)))
                d = {"hazard":hz,"sector":s}
                add_sectors.append(dict(d,**val_dict))
        add_sectors = pd.DataFrame(add_sectors)
        df_hz = pd.concat([df_hz,add_sectors],axis=0,ignore_index=True)
        df_hz = df_hz.sort_values(by=["sector"],ascending=True)

        # print (hz)
        print (df_hz)
        fig, ax = plt.subplots(1,1,figsize=(8,8),dpi=500)
        ax.barh(df_hz["sector"],multiply_factor*df_hz["EAD_mean"],color=d_color,label="Direct damages (mean value)")
        ax.barh(df_hz["sector"],multiply_factor*df_hz["EAEL_mean"],color=i_color,left=multiply_factor*df_hz["EAD_mean"],
                    xerr=(multiply_factor*(df_hz["EAD_mean"] - df_hz["EAD_amin"]),multiply_factor*(df_hz["EAD_amax"] - df_hz["EAD_mean"])),
                    capsize=5,label="Indirect economic losses (mean value)")
        
        # ax.legend(prop={'size':12,'weight':'bold'})
        ax.legend(handles=legend_handles,prop={'size':12,'weight':'bold'})
        ax.set_xlabel('Expected Annual Damages and Losses (J$ Billion)',fontweight='bold',fontsize=15)
        ax.set_xlim(0,max_value)
        # ax.set_xscale('log')
        ax.set_yticklabels(df_hz["sector"].values.tolist(),fontsize=15, rotation=0,fontweight="bold")
        ax.set_title(f"{hz_l} Direct and Indirect Risks",fontweight="bold",fontsize=15)
        plt.tight_layout()
        save_fig(os.path.join(figures_data_path,
                    f"jamaica_hazard_EAD_EAEL_{hz}_totals.png"))
        plt.close()


    legend_handles = []
    legend_handles.append(mpatches.Patch(color=d_color, label="Direct Damages (Mean Value)"))
    legend_handles.append(Line2D([0], [0], 
                            color=l_color, label="Min-Max Range"))

    df = ead_results[ead_results["epoch"] == 2010].groupby(["hazard","sector"])[ead_columns].sum().reset_index()
    all_sectors = list(set(df["sector"].values.tolist()))
    print (df)
    max_value = multiply_factor*df["EAD_amax"].max() + 2

    for idx,(hz,hz_l) in enumerate(list(zip(hazards,hazard_labels))):
        df_hz = df[df["hazard"] == hz]
        add_sectors = []
        for s in all_sectors:
            if s not in df_hz["sector"].values.tolist():
                val_dict = dict(zip(ead_columns,[0]*len(ead_columns)))
                d = {"hazard":hz,"sector":s}
                add_sectors.append(dict(d,**val_dict))
        add_sectors = pd.DataFrame(add_sectors)
        df_hz = pd.concat([df_hz,add_sectors],axis=0,ignore_index=True)
        df_hz = df_hz.sort_values(by=["sector"],ascending=True)

        fig, ax = plt.subplots(1,1,figsize=(8,8),dpi=500)
        ax.barh(df_hz["sector"],multiply_factor*df_hz["EAD_mean"],color=d_color,
                    xerr=(multiply_factor*(df_hz["EAD_mean"] - df_hz["EAD_amin"]),multiply_factor*(df_hz["EAD_amax"] - df_hz["EAD_mean"])),
                    capsize=5,label="Direct damages (mean value)")
        
        # ax.legend(prop={'size':12,'weight':'bold'})

        ax.legend(handles=legend_handles,prop={'size':12,'weight':'bold'})
        ax.set_xlabel('Expected Annual Damages (J$ Billion)',fontweight='bold',fontsize=15)
        ax.set_xlim(0,max_value)
        ax.set_yticklabels(df_hz["sector"].values.tolist(),fontsize=15, rotation=0,fontweight="bold")
        ax.set_title(f"{hz_l} Direct Physical Risks",fontweight="bold",fontsize=15)
        plt.tight_layout()
        save_fig(os.path.join(figures_data_path,
                    f"jamaica_hazard_EAD_{hz}_totals.png"))
        plt.close()


    """Timeseries plots
    """
    timeseries = pd.read_csv(os.path.join(damage_data_path,
                                    "hazard_rcp_timeseries_totals.csv"))
    # timeseries = timeseries[(timeseries["sector"] != "Potable water") | (timeseries["subsector"] != "edges")]
    start_year = 2019
    end_year = 2100
    years = np.arange(start_year,end_year+1,1)
    damage_columns = [str(n) for n in np.arange(start_year,end_year+1,1)]

    # index_columns = ["hazard","rcp","risk_type","val_type"]
    index_columns = ["hazard","rcp","val_type"]
    df = timeseries.groupby(index_columns)[damage_columns].sum().reset_index()
    max_value = multiply_factor*max(df[damage_columns].max())

    rcp_colors = ["#31a354","#3182bd","#de2d26"]
    rcps = [2.6,4.5,8.5]
    rcp_markers= ['p-','*-','o-']
    for ix,(hz,hz_l) in enumerate(list(zip(hazards,hazard_labels))):
        fig, ax = plt.subplots(1,1,figsize=(8,8),dpi=500)
        for idx, (rcp,rcp_c,rcp_m) in enumerate(list(zip(rcps,rcp_colors,rcp_markers))):
            rcp_df = df[(df["hazard"] == hz) & (df["rcp"] == rcp)]
            if len(rcp_df.index) > 0:
                mean_vals = rcp_df[rcp_df["val_type"] == "mean"][damage_columns].values[0]
                min_vals = rcp_df[rcp_df["val_type"] == "amin"][damage_columns].values[0]
                max_vals = rcp_df[rcp_df["val_type"] == "amax"][damage_columns].values[0]
                # print (mean_vals)
                ax.plot(years,multiply_factor*mean_vals,rcp_m,color=rcp_c,markersize=4,linewidth=2.0,label=f"RCP {rcp} mean")
                ax.fill_between(years,multiply_factor*min_vals,multiply_factor*max_vals,alpha=0.3,facecolor=rcp_c,label=f"RCP {rcp} min-max")

        ax.legend(prop={'size':12,'weight':'bold'},loc="upper left")
        ax.set_ylabel('Expected Annual Damages + Losses (J$ Billion)',fontweight='bold',fontsize=15)
        ax.set_xlabel('Years',fontweight='bold',fontsize=15)
        ax.set_ylim(0,max_value)
        ax.set_title(f"{hz_l} Direct and Indirect Risks over time",fontweight="bold",fontsize=15)
        plt.tight_layout()
        save_fig(os.path.join(figures_data_path,
                    f"jamaica_hazard_risk_{hz}_timeseries.png"))
        plt.close()


if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
