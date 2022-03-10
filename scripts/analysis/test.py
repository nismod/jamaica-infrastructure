"""Functions for preprocessing data
"""
import sys
import os
import json

import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import geopandas as gpd
from scipy.interpolate import interp1d
from scipy import integrate
import math
import numpy as np
from tqdm import tqdm
tqdm.pandas()
from analysis_utils import *


def main(config):
    A = np.array([1,2,3,4])
    B = np.array([1,4,6,8])
    df = pd.DataFrame()
    df["W"] = [-1,1.5,3,2.5,4.5,3.5]
    df["X"] = [0,1,2,4,5,6]
    print (df)
    df[["Y","Z"]] = interp1d(A,B,fill_value=(min(B),max(B)),bounds_error=False)(df[["W","X"]])
    print (df)
    df[["Y","Z"]] = df[["Y","Z"]] - 1
    df[["Y","Z"]] = df["X"].to_numpy()[:,None]*np.where(df[["Y","Z"]] > 0,1,0)
    print (df)
    df = df[(df[["Y","Z"]]>4).all(axis=1)]
    print (df)
    # df[df[["Y","Z"]] < 7] = 10
    # print (df)
    df[["Y","Z"]] = np.where(df[["Y","Z"]]<7,df[["Y","Z"]],10)
    print (df)

    A = np.array([[0.01,0.1]*len(df.index)])
    print (A)
    df["A"] = integrate.trapz(df[["Y","Z"]].to_numpy(),A)
    print (df)
    df[["Y","Z"]] = ((df["X"].to_numpy().T)*np.where(df[["Y","Z"]] > 0,1,0)).T
    print (df)
if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)