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
mpl.rcParams['xtick.labelsize'] = 10.
mpl.rcParams['ytick.labelsize'] = 10.


def main(config):
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    figures_data_path = config['paths']['figures']

    damage_data_path = os.path.join(incoming_data_path,
                        'damage_assessment_information',
                        'damage_functions')

    damage_data_keys = pd.read_csv(os.path.join(damage_data_path,'damage_curves_keys.csv'))

    for data_key in damage_data_keys.itertuples():
        data = pd.read_excel(os.path.join(damage_data_path,
                            f"damage_curves_{data_key.sector}_{data_key.hazard_type}.xlsx"),
                            sheet_name=data_key.asset_sheet)

        fig, ax = plt.subplots(1,1,figsize=(8,8),dpi=500)
        if data_key.hazard_type == 'flooding':
            x_data = data['flood_depth']
            x_label = 'Flood depth (m)'
        else:
            x_data = data['wind_speed']
            x_label = 'Wind speed (m/s)'

        if 'damage_ratio_1' in data.columns.values.tolist():
            ax.plot(x_data,
                    data['damage_ratio_1'],'-',
                    color='#000000',linewidth=2.0,
                    label=f'{data_key.source_1} curve')
            ax.plot(x_data,
                    data['damage_ratio_2'],'-',
                    color='#cb181d',linewidth=2.0,
                    label=f'{data_key.source_2} curve')
        else:
            ax.plot(x_data,
                    data['damage_ratio'],'-',
                    color='#3182bd',linewidth=2.0,
                    label=f'{data_key.source_1} curve')    

        ax.fill_between(x_data,
                        data['damage_ratio_min'],
                        data['damage_ratio_max'],
                        alpha=0.3,facecolor='#3182bd',
                        label='Uncertainty range')
        ax.legend(loc='upper left',fontsize=14)
        ax.set_xlabel(x_label,fontweight='bold',fontsize=12)
        ax.set_ylabel('Damage ratio [0-1]',fontweight='bold',fontsize=12)
        if data_key.hazard_type == 'TC':
           title_label = 'Wind damage'
        else:
            title_label = 'Flood damage'

        ax.set_title(f"{data_key.asset_label}: {title_label}",fontsize=14,fontweight='bold')
        plt.tight_layout()
        save_fig(os.path.join(figures_data_path,
                    f'damage_curves-{data_key.hazard_type}-{data_key.asset_sheet}.png'))
        plt.close()    

if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
