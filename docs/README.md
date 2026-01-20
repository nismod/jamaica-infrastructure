# A geospatial analysis platform for infrastructure risk assessment and resilient investment prioritisation in Jamaica

## Methodology and implementation of the Jamaica Systemic Risk Assessment Tool (J-SRAT)

This document describes the methodology, data and implementation for the
development of the Jamaica Systemic Risk Assessment Tool (J-SRAT). The
focus of this report is to:

1.  Provide a background of the project;

2.  Provide a detailed overview of the methodology and the component
    models;

3.  Outline the modelling assumptions made in the analysis;

4.  Describe the datasets collected and finalised in implementing the
    methodology;

5.  Show the process of implementation of the methodology through a use
    case example.

This is a _technical document_ which is meant to explain the theory and
concepts behind the spatial risk analysis and component models being
developed in Jamaica. The _main target audience_ of this report include
climate risk and infrastructure modellers and technical experts who
would be interested in understanding the detailed working of the J-SRAT.

## Contents

1. [Introduction](01-introduction.md)

   - Background
   - The need and opportunity for infrastructure risk and resilience assessment
   - Project objectives and impact

2. [Methodology development and implementation steps](02-methodology.md)

   - Types of assessment done through J-SRAT
   - J-SRAT framework details and risk calculations
   - Output metrics
   - Implementation steps of the J-SRAT methodology

3. [Model assumptions and data assembled for implementing J-SRAT](03-model-assumptions-data.md)

   - Generic methodology assumptions
   - Hazard data assembly
   - Socio-economic model and data assembly
   - Infrastructure asset and network model and data assembly
   - Adaptation options

## Appendices

- [Appendix A: Vulnerability curves for infrastructure assets in Jamaica](appendix-a-vulnerability-curves.md)

- [Appendix B: Hazard models](appendix-b-hazard-models.md)

  - B.1 Fluvial and pluvial flood models
  - B.2 Coastal flood model
  - B.3 Tropical cyclone model
  - B.4 Drought model

- [Appendix C: Infrastructure network flow models for failure analysis](appendix-c-network-flow-models.md)

  - C.1 Generalised network representations
  - C.2 Energy system model
  - C.3 Water systems models
  - C.4 Transport system model

- [Appendix D: Spatial disaggregation of economic activity at buildings and area levels](appendix-d-spatial-disaggregation.md)
  - D.1 National accounting of GDP by economic sectors and subsectors
  - D.2 Disaggregation of GDP to buildings
  - D.3 Disaggregating GDP to agriculture areas
  - D.4 Disaggregating GDP to mining and quarrying areas

## Document control

| Issue | Status                                  | Author(s)        | Reviewed by | Issue Date |
| ----- | --------------------------------------- | ---------------- | ----------- | ---------- |
| 1     | Final report                            | Raghav Pant      | Tim Fowler  | Jul 2022   |
|       |                                         | Olivia Becher    | Jim Hall    |            |
|       |                                         | Robyn Haggis     |             |            |
|       |                                         | Aman Majid       |             |            |
|       |                                         | Tom Russell      |             |            |
|       |                                         | Jasper Verschuur |             |            |
|       |                                         | Edson Williams   |             |            |
|       |                                         | Anaitée Mills    |             |            |
|       |                                         | Ardith Grant     |             |            |
|       |                                         |                  |             |            |
| 2     | Conversion to markdown documentation    | as above, with:  |             | Jan 2026   |
|       | and focus on methods, not data/results. | Fred Thomas      |             |            |

This report may be cited as follows:

> Pant, R., Becher, O., Haggis R., Majid, A., Russell, T., Verschuur, J.,
> Williams, E., Mills, A., Grant, A. Fowler, T. and Hall, J.W. (2022). Final
> technical report on methodology and implementation of the Jamaica Systemic
> Risk Assessment Tool (J-SRAT). Environmental Change Institute, Oxford
> University, UK.

Initial Report © Oxford University, 2022

This report was produced as part of the project "A geospatial analysis platform
for infrastructure risk assessment and resilient investment prioritisation in
Jamaica" which was funded by the UK Foreign Commonwealth and Development Office
(FDCO) as part of the Coalition for Climate Resilient Investment (CCRI). The
views expressed and recommendations set out in this report are the authors' own
and do not necessarily reflect the position of the FCDO, CCRI or any public or
private stakeholder in Jamaica.

The materials have been prepared by Oxford University. Whilst every care has
been taken by Oxford University to ensure the accuracy and completeness of the
reports and maps, the reader must recognise that errors are possible through no
fault of Oxford University and as such the parties give no express or implied
representations or warranty as to:

\(i\) the quality or fitness for any particular purpose of the report or
maps supplied or of any design, workmanship, materials or parts used in
connection therewith or correspondence with regard to any description or
sample; or

\(ii\) the accuracy, sufficiency or completeness of the reports or maps
provided. In particular, there are hereby expressly excluded all
conditions, warranties and other terms which might otherwise be implied
(whether by common law, by statute or otherwise).

Oxford University, its employees, servants and agents shall accept no
liability for any damage caused directly or indirectly by the use of any
information contained herein and without prejudice to the generality of
the foregoing, by any inaccuracies, defects or omissions.

![Logos of OPSIS, Environmental Change Institute (ECI), and University of Oxford](media/oxford-eci-opsis-logos.png)
