"""Create population projections for Jamaica
"""
import sys
import os

import pandas as pd
import geopandas as gpd
import numpy as np

from utils import *
from tqdm import tqdm
tqdm.pandas()

def get_age_bins(x):
    if x["age_group"] == "65 and over":
        return 65,69
    else:
        return int(x["age_group"].split("-")[0]),int(x["age_group"].split("-")[1])

def modify_working_ages(working_dataframe):
    df = []
    for row in working_dataframe.itertuples():
        if row.min_age == 14:
            df.append(("15-19",row.min,row.mean,row.max))
        elif row.max_age - row.min_age == 4:
            df.append((f"{row.min_age}-{row.max_age}",row.min,row.mean,row.max))
        else:
            df.append((f"{row.min_age}-{row.min_age+4}",row.min,row.mean,row.max))
            df.append((f"{row.max_age-4}-{row.max_age}",row.min,row.mean,row.max))

    return pd.DataFrame(df,columns=["age_group","percent_employed_min","percent_employed_mean","percent_employed_max"]) 

def modify_population_age_change(population_dataframe):
    df = []
    for row in population_dataframe.itertuples():
        if "-" not in row.age_group:
            df.append(("80-84",row.growth_rate_min,row.growth_rate_mean,row.growth_rate_max))
            df.append(("85-89",row.growth_rate_min,row.growth_rate_mean,row.growth_rate_max))
            df.append(("90-94",row.growth_rate_min,row.growth_rate_mean,row.growth_rate_max))
            df.append(("95-",row.growth_rate_min,row.growth_rate_mean,row.growth_rate_max))
        else:
            df.append((row.age_group.strip(),row.growth_rate_min,row.growth_rate_mean,row.growth_rate_max))

    return pd.DataFrame(df,columns=["age_group","growth_rate_min","growth_rate_mean","growth_rate_max"]) 

def get_population_changes(population_dataframe,start_year,end_year,base_year):
    growth_rate_columns = []
    for year in range(start_year,end_year):
        year_column = str(year)
        base_year_column = str(base_year)
        population_dataframe[year_column] = population_dataframe.progress_apply(lambda x:int(str(x[year_column]).replace(",","")),axis=1)
        population_dataframe[f"growth_rate_{year}"] = (np.power(population_dataframe[year_column]/population_dataframe[base_year_column],1/(year-base_year))-1)
        growth_rate_columns.append(f"growth_rate_{year}")

    population_dataframe["growth_rate_mean"] = population_dataframe[growth_rate_columns].mean(axis=1)
    population_dataframe["growth_rate_min"] = population_dataframe[growth_rate_columns].min(axis=1)
    population_dataframe["growth_rate_max"] = population_dataframe[growth_rate_columns].max(axis=1)
    return population_dataframe

def main(config):
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    epsg_jamaica = 3448

    baseyear = 2011
    id_column = "ED_ID"
    id_columns = ["ED_ID","ED"]
    other_columns = ['ED_CLASS','AREA', 'PERIMETER', 'PARISH', 'CONST_NAME', 'ED']
    male_age_wise_columns = [
        'F0_4_MLE', 'F5_9_MLE', 'F10_14_MLE', 'F15_19_MLE', 'F20_24_MLE',
        'F25_29_MLE', 'F30_34_MLE', 'F35_39_MLE', 'F40_44_MLE', 'F45_49_MLE',
        'F50_54_MLE', 'F55_59_MLE', 'F60_64_MLE', 'F65_69_MLE', 'F70_74_MLE',
        'F75_79_MLE', 'F80_84_MLE', 'F85_89_MLE', 'F90_94_MLE', 'F95__MLE']
    male_working_columns = [
        'F15_19_MLE', 'F20_24_MLE',
        'F25_29_MLE', 'F30_34_MLE', 'F35_39_MLE', 'F40_44_MLE', 'F45_49_MLE',
        'F50_54_MLE', 'F55_59_MLE', 'F60_64_MLE', 'F65_69_MLE','F70_74_MLE',
        'F75_79_MLE']
    male_nonworking_columns = [c for c in male_age_wise_columns if c not in male_working_columns]
    female_age_wise_columns = [
        'F0_4_FMLE', 'F5_9_FMLE', 'F10_14_FMLE', 'F15_19_FMLE',
        'F20_24_FMLE', 'F25_29_FMLE', 'F30_34_FMLE', 'F35_39_FMLE',
        'F40_44_FMLE', 'F45_49_FMLE', 'F50_54_FMLE', 'F55_59_FMLE',
        'F60_64_FMLE', 'F65_69_FMLE', 'F70_74_FMLE', 'F75_79_FMLE',
        'F80_84_FMLE', 'F85_89_FMLE', 'F90_94_FMLE', 'F95__FMLE']
    female_working_columns = [
        'F15_19_FMLE',
        'F20_24_FMLE', 'F25_29_FMLE', 'F30_34_FMLE', 'F35_39_FMLE',
        'F40_44_FMLE', 'F45_49_FMLE', 'F50_54_FMLE', 'F55_59_FMLE',
        'F60_64_FMLE', 'F65_69_FMLE','F70_74_FMLE', 'F75_79_FMLE']
    female_nonworking_columns = [c for c in female_age_wise_columns if c not in male_working_columns]

    female_total_column = "TOTAL_FMLE"
    male_total_column = "TOTAL_MLE"
    total_column_column = 'TOTAL_POP' 
    population = gpd.read_file(os.path.join(processed_data_path,
                            'boundaries',
                            'admin_boundaries.gpkg'),layer='admin3')
    population_estimation = gpd.read_file(os.path.join(processed_data_path,
                            'population',
                            'population.gpkg'),layer='admin3')
    population_estimation["parish_rename"] = population_estimation.progress_apply(lambda x:str(x["PARISH"]).replace(".","").replace(" ","").lower(),axis=1)

    population_growth = pd.read_excel(os.path.join(incoming_data_path,
                                    "macroeconomic_data",
                                    "parish_population_changes.xlsx"),
                                    sheet_name="2014-2019")
    population_growth = population_growth[population_growth["Parish"] != "Total"]
    population_growth["parish_rename"] = population_growth.progress_apply(lambda x:str(x["Parish"]).replace(" ","").lower(),axis=1)

    population = pd.merge(population,population_estimation[id_columns+["population"]],how="left",on=id_columns)
    population["parish_rename"] = population.progress_apply(lambda x:str(x["PARISH"]).replace(".","").replace(" ","").lower(),axis=1)
    """Get the population growth from 2011 to 2019
    """
    parish_population = population.groupby("parish_rename")["population"].sum().reset_index()
    parish_population_change = pd.merge(parish_population,population_growth,how="left",on=["parish_rename"])
    parish_population_change.rename(columns={"population":"2011"},inplace=True)
    parish_population_change.columns = parish_population_change.columns.map(str)
    parish_population_change = get_population_changes(parish_population_change,2014,2020,baseyear)
    print (parish_population_change)
    
    """Population changes by age group
    """
    forecasts = ["min","mean","max"]
    employed_population = pd.read_csv(os.path.join(incoming_data_path,
                                    "macroeconomic_data",
                                    "population_employment_percent_by_age.csv"))
    employed_population["age_group"] = employed_population.apply(lambda x:str(x["age_group"]).strip().replace(" - ","-"),axis=1)
    employed_population = employed_population[employed_population["age_group"] != "TOTAL"]
    working_categories = [("Male","male_age_groups.csv","MLE"),("Female","female_age_groups.csv","FMLE")]
    pop_totals = []
    all_year_columns = []
    for i,(gender,gender_file,gender_string) in enumerate(working_categories):
        gender_working = employed_population[employed_population["catetorgy"] == gender]
        gender_working["age_bins"] = gender_working.progress_apply(lambda x:get_age_bins(x),axis=1)
        gender_working[["min_age","max_age"]] = gender_working["age_bins"].apply(pd.Series)
        population_changes = modify_working_ages(gender_working)
        population_gender = population.copy()
        for forecast in forecasts:
            year_column = f"{forecast}_{gender.lower()}_{baseyear}"
            all_year_columns.append(year_column)
            population_gender[year_column] = 0
            for pop in population_changes.itertuples():
                pop_column = f"F{str(pop.age_group).replace('-','_')}_{gender_string}"
                population_gender[year_column] += population_gender[pop_column]*0.01*getattr(pop,f"percent_employed_{forecast}")

            population_gender[year_column] = population_gender.apply(lambda x:round(x[str(year_column)],0),axis=1)

        population_gender.drop(male_age_wise_columns+female_age_wise_columns+["parish_rename","population"],axis=1,inplace=True)
        pop_totals.append(population_gender[[c for c in population_gender.columns.values.tolist() if c != "geometry"]])

    index_cols = [c for c in pop_totals[0].columns.values.tolist() if c not in all_year_columns]
    pop_totals = pd.merge(pop_totals[0],pop_totals[1],how="left",on=index_cols)
    for forecast in forecasts:
        pop_totals[f"{forecast}_working_{baseyear}"] = pop_totals[f"{forecast}_male_{baseyear}"] + pop_totals[f"{forecast}_female_{baseyear}"] 

    print (pop_totals)
    print (population_estimation)
    """Project population to the end of the century
    """
    population_estimation = pd.merge(population_estimation,
                            parish_population_change[["parish_rename","growth_rate_min",
                                                    "growth_rate_mean","growth_rate_max"]],
                            how="left",on=["parish_rename"])
    population_estimation.rename(columns={"population":f"{baseyear}"},inplace=True)
    population_estimation = pd.merge(population_estimation,
                                pop_totals[id_columns + [f"min_working_{baseyear}",
                                                f"mean_working_{baseyear}",
                                                f"max_working_{baseyear}"]],
                                how="left",on=id_columns)
    population_estimation.drop("parish_rename",axis=1,inplace=True)
    forecasts = ["min","mean","max"]
    exclude_columns = ["growth_rate_min",
                        "growth_rate_mean",
                        "growth_rate_max",
                        f"min_working_{baseyear}",
                        f"mean_working_{baseyear}",
                        f"max_working_{baseyear}"]
    for forecast in forecasts:
        population_estimation["working_frac"] = population_estimation.apply(
                                        lambda x:x[f"{forecast}_working_{baseyear}"]/x["2011"] if x["2011"] > 0 else 0,
                                        axis=1)
        if f"working_{baseyear}" in population_estimation.columns.values.tolist():
            population_estimation.drop(f"working_{baseyear}",axis=1,inplace=True)

        population_estimation.rename(columns={f"{forecast}_working_{baseyear}":f"working_{baseyear}"},inplace=True)
        for year in range(2012,2100):
            population_estimation[str(year)] = population_estimation[str(baseyear)]*np.power(
                                        1+population_estimation[f"growth_rate_{forecast}"],
                                        year-baseyear)
            population_estimation[f"working_{year}"] = population_estimation["working_frac"]*population_estimation[str(year)]
            population_estimation[str(year)] = population_estimation.apply(lambda x:round(x[str(year)],0),axis=1)
            population_estimation[f"working_{year}"] = population_estimation.apply(lambda x:round(x[f"working_{year}"],0),axis=1)

        population_estimation.drop("working_frac",axis=1,inplace=True)
        result_gdf = population_estimation[[c for c in population_estimation.columns.values.tolist() if c not in exclude_columns]]
        result_gdf.to_file(os.path.join(processed_data_path,
                            'population',
                            'population_projections.gpkg'),layer=forecast,driver="GPKG")

if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)