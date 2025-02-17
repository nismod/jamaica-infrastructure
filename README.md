# Jamaica Infrastructure Risk and Resilience Assessment

This repository contains project-specific codes and configuration to run
climate-related risk and resilience analysis of infrastructure networks in
Jamaica.

For an initial overview of the style of analysis, a table of potential data
requirements is presented in [data-categories.csv](./data-categories.csv).

## Setup and installation

Clone or download this repository from
[GitHub](https://github.com/nismod/jamaica-infrastructure):

    git clone git@github.com:nismod/jamaica-infrastructure.git

Next, install required python packages:

We recommend using [`micromamba`](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)
to install the packages and manage installing libraries into a conda
environment, usually handling non-Python dependencies well.

Create a conda environment once (per machine/user):

    micromamba create --file environment.yml

Activate the environment each time you open a shell:

    micormamba activate jsrat

## High-level analysis steps

Run hazard-network intersections.

```bash
python scripts/exposure/split_networks.py \
    workflow/network_layers.csv \
    workflow/hazard_layers.csv \
    ./extract/processed_data/
```



## Related repositories

### Energy

The [Jamaica Energy Model (JEM)](https://github.com/nismod/jem) is a high-level
power flow model of Jamaica's electricity network.

### Visualisation tool

The J-SRAT web-based visualisation tool is implemented in
[`irv-jamaica`](https://github.com/nismod/irv-jamaica)

### snail

[`snail`](https://github.com/nismod/snail) is a supporting library, used here
for hazard-network intersections, under continuing development to support risk
analysis.

### open-gira

[`open-gira`](https://github.com/nismod/open-gira) is a supporting library and
set of workflows, used here for road network generation

## Acknowledgments

This work is supported by the
[Coalition for Climate Resilient Investment (CCRI)](https://resilientinvestment.org/)
project on infrastructure risk assessment and resilient investment
prioritisation in Jamaica, funded by the
[UK Foreign, Commonwealth and Development Office (FCDO)](https://www.gov.uk/government/organisations/foreign-commonwealth-development-office).
