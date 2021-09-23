#!/usr/bin/env python
# coding: utf-8
"""Parse hazard layer metadata from path/filenames
"""
import re

import pandas

# Read basic info (createcd from running
# `ls hazards/**/**/**/*.tif >> hazard_layers_basic.csv` and annotating with
# hazard type)
df = pandas.read_csv("hazard_layers_basic.csv")

# Set return period
df["rp"] = None

df.loc[df.hazard.isin(("fluvial", "surface")), "rp"] = df[
    df.hazard.isin(("fluvial", "surface"))
].path.apply(lambda p: int(re.search(r"Q(\d+)", p).group(1)))

df.loc[df.hazard == "coastal", "rp"] = df[df.hazard == "coastal"].path.apply(
    lambda p: int(re.search(r"RP_(\d+)", p).group(1))
)

df.loc[df.hazard == "cyclone", "rp"] = df[df.hazard == "cyclone"].path.apply(
    lambda p: int(re.search(r"rp(\d+)", p).group(1))
)

# Set RCP
df["rcp"] = None

df.loc[df.hazard.isin(("fluvial", "surface")), "rcp"] = "baseline"

df.loc[df.hazard == "coastal", "rcp"] = df[df.hazard == "coastal"].path.apply(
    lambda p: int(re.search(r"RCP(\d\d)", p).group(1)) / 10
)


def cyclone_rcp(p):
    if "baseline" in p:
        return "baseline"
    return int(re.search(r"RCP(\d\d)", p).group(1)) / 10


df.loc[df.hazard == "cyclone", "rcp"] = df[df.hazard == "cyclone"].path.apply(
    cyclone_rcp
)

# Set epoch
df["epoch"] = None

df.loc[df.hazard.isin(("fluvial", "surface")), "epoch"] = 2010

df.loc[df.hazard == "coastal", "epoch"] = df[df.hazard == "coastal"].path.apply(
    lambda p: int(re.search(r"RCP\d\d(\d\d\d\d)", p).group(1))
)


def cyclone_epoch(p):
    if "baseline" in p:
        return 2010
    if "midcentury" in p:
        return 2050
    if "endcentury" in p:
        return 2100


df.loc[df.hazard == "cyclone", "epoch"] = df[df.hazard == "cyclone"].path.apply(
    cyclone_epoch
)

# Set confidence
df["confidence"] = None


def cyclone_confidence(p):
    if "_mean" in p:
        return 50
    if "_conf_5" in p:
        return 5
    if "_conf_95" in p:
        return 95


df.loc[df.hazard == "cyclone", "confidence"] = df[df.hazard == "cyclone"].path.apply(
    cyclone_confidence
)

# Sort
df = df.sort_values(by=["hazard", "rcp", "epoch", "rp", "confidence"]).reset_index(
    drop=True
)

# Set key (encode all values)
df["key"] = df.apply(
    lambda h: f"{h.hazard}__rp_{h.rp}__rcp_{h.rcp}__epoch_{h.epoch}__conf_{h.confidence}",
    axis=1,
)

# Save
df.to_csv("hazard_layers.csv")
