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
```shell
micromamba create --file environment.yml
```

## Usage

The principal goal of this repository is to produce analysis results, (damages,
losses and adaptation options) that can be ingested by the ETL pipeline in
[irv-jamaica](https://github.com/nismod/irv-jamaica/blob/main/etl/README.md).

The analysis is comprised of Python scripts wrapped in
[snakemake](https://snakemake.readthedocs.io/) rules. Snakemake is a workflow
management system that can be used to break up complex modelling chains and
improve the reproducibility of analyses. To invoke a rule, call `snakemake`
followed by the file you want to produce.

First, make available snakemake and other software dependencies by activating
the environment we previously created.
```shell
micromamba activate jsrat
```

And to, for example, invoke the rule (and all necessary predecessor rules) to
compute commuter flows across the transport network:
```shell
snakemake --dry-run --cores 1 -- results/flow_mapping/labour_to_sectors_flow_paths.csv
```

Note that the `--dry-run` flag asks `snakemake` to report on what work (if any)
is necessary. As jobs can be long running, this is useful to check beforehand.
To actually run the rules, remove the dry run flag.

The `--cores` flag indicates how many processors `snakemake` will use to execute
the rules. If rules do not depend on one another and enough processors are
available, they may execute simultaneously. Also, some rules invoke scripts that
are parallelised and can make use more than one processor themselves.

### Rules

See the following files within `workflow/` for available rules and their input
and output files:
- `direct_damages.smk`
- `losses.smk`
- `transport_model.smk`
- `hotspots.smk`

## Related repositories

### Energy

The [Jamaica Energy Model (JEM)](https://github.com/nismod/jem) is a high-level
power flow model of Jamaica's electricity network.

### Visualisation tool

The J-SRAT web-based visualisation tool is implemented in
[`infra-risk-vis`](https://github.com/nismod/infra-risk-vis)

### snail

[`snail`](https://github.com/nismod/snail) is a supporting library, used here
for hazard-network intersections, under continuing development to support risk
analysis.

## Acknowledgments

This work is supported by the
[Coalition for Climate Resilient Investment (CCRI)](https://resilientinvestment.org/)
project on infrastructure risk assessment and resilient investment
prioritisation in Jamaica, funded by the
[UK Foreign, Commonwealth and Development Office (FCDO)](https://www.gov.uk/government/organisations/foreign-commonwealth-development-office).
