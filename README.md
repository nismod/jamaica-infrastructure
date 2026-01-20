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

Next, install required Python packages:

We recommend using [`micromamba`](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)
to install the packages and manage installing libraries into a conda
environment, usually handling non-Python dependencies well.

Create a conda environment once (per machine/user):

    micromamba create --file environment.yaml

## Usage

The principal goal of this repository is to produce analysis results, (damages,
losses and adaptation options) that can be ingested by the ETL pipeline in
[irv-jamaica](https://github.com/nismod/irv-jamaica/blob/main/etl/README.md).

The analysis is comprised of Python scripts wrapped in
[snakemake](https://snakemake.readthedocs.io/) rules. Snakemake is a workflow
management system that can be used to break up complex modelling chains and
improve the reproducibility of analyses. To invoke a rule, call `snakemake`
followed by the file you want to produce.

Some functionality is contained within a helper Python module, located in
`src/jamaica_infrastructure`.

### Workflow overview

The analysis pipeline is implemented using Snakemake and organized into several main stages:

1. **Damage Assessment** (`1_damage`) - Calculate direct damages to infrastructure assets from hazards
2. **Transport Analysis** (`2a_transport`) - Create multi-modal transport networks and flow mapping
3. **Criticality Analysis** (`3_criticality`) - Assess network performance under single-point failures
4. **Losses** (`4a_losses`) - Calculate Expected Annual Damages (EAD) and Expected Annual Economic Losses (EAEL)
5. **Hotspots** (`4b_hotspots`) - Identify spatial concentrations of risk
6. **Self-Adaptation** (`5a_self-adaptation`) - Evaluate asset-level adaptation options and benefit-cost ratios
7. **Coastal Adaptation** (`5b_coastal_adaptation`) - Generate coastal protection scenarios
8. **ETL** (`6_etl`) - Prepare data for visualization

### Environment

To make snakemake, the helper Python module and other software dependencies
available, activate the environment we previously created.
```shell
micromamba activate jsrat
```

### Configuration

Analysis options can be configured using the `config.yaml` file. See it for
inline documentation.

### Invoking rules

To invoke the rule (and all necessary predecessor rules) to compute commuter
flows across the transport network:
```shell
snakemake --dry-run --cores 1 -- results/flow_mapping/labour_to_sectors_flow_paths.pq
```

Note that the `--dry-run` flag asks `snakemake` to report on what work (if any)
is necessary. As jobs can be long running, this is useful to check beforehand.
To actually run the rules, remove the dry run flag.

The `--cores` flag indicates how many processors `snakemake` will use to execute
the rules. If rules do not depend on one another and enough processors are
available, they may execute simultaneously. Also, some rules invoke scripts that
are parallelised and can make use more than one processor themselves.

### Tests

The helper Python libary contained in `src/jamaica_infrastructure` also has
tests. These can be run with:

    python -m pytest src/jamaica_infrastructure

These take a few seconds.

### Required data

To run any rules you will need the `processed_data/` folder which contains data
that has been cleaned and can be consumed by the rules. Contact the maintainers
for access to this folder.

Most(!) outputs are written to the `results/` folder.

These paths can be configured by editing the `config.yaml` file.

More information on the input and output data can be found in `DATA.md`.

### Development

Outstanding work may be seen the [issues](https://github.com/nismod/jamaica-infrastructure/)
on the GitHub repository.

If you are resurrecting scripts pertaining to this repository, they may rely on
functionality since removed. Checkout v1.0 for the state of the project as of
2023 (prior to the rewrite of this repository in phase 3).

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
