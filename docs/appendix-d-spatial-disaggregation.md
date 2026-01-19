# Appendix D: Spatial disaggregation of economic activity at buildings and area levels

The spatial disaggregation of the economic activity in Jamaica involved
taking national-scale statistics and creating methods to disaggregate
them to individual buildings and areas. The steps involved in the
process are described in the following sections.

## D.1 National accounting of GDP by economic sectors and subsectors

1.  We first obtained the national accounts of sector specific Gross
    Value Added (GVA) and total GDP for 2019 as reported by STATIN
    (<https://statinja.gov.jm/NationalAccounting/Annual/NewAnnualGDP.aspx>,
    see Table D-1).

2.  From the STATIN data we can see the GVA + TAX = GDP, and hence we
    can estimate the percentage of TAX added to GVA to get GDP. From
    Table D-1

TAX = 23% of GVA (D.1)

3.  We assume that the TAX contribution to each sector to the economic
    is in the same percentage, which implies that for a given sector
    _s_:

(1 + 0.23)×GVA~s~ = GDP~s~ (D.2)

## D.2 Disaggregation of GDP to buildings

1.  For the sectors D -- O, listed in Table 3‑5 and Table D-1, we
    created the geotagged buildings database described in Section 3.3.2.
    We assume that the total sector GDP will be disaggregated to these
    buildings.

2.  The amount of GDP associated with a building is estimated to be a
    function of the location of the buildings in proximity with working
    populations and the area of the building. The proximity to working
    population is a measure of the attractiveness of the area where the
    building is located and the area is a measure of the capacity of the
    building in accommodating more workforce.

3.  To estimate the attractiveness of the area where a building is
    located we used the Enumeration District (ED) level total and
    working population estimates, from the data described in Section
    3.3.1. We then apply the radiation model (introduced in Section C.4)
    to infer how mobility patterns will determine concentrations of GDP.
    For the collection of ED areas in Jamaica we know the total
    populations $p_{i},p_{j}$ and employed populations $e_{i},e_{j}$
    between two ED-pairs $(i,j)$ located within a distance (or travel
    time) of $d_{ij}$. The number of people $T_{ij}$ who will move for
    ED-area $i$ to ED-area $j$ to seek employment is estimated from a
    radiation model^76^ as shown in Equation D.3. Here $e_{ij}$ is equal
    to the employed population within the radius defined by the travel
    distance $d_{ij}$.

$T_{ij} = \ p_{i}\frac{e_{j}/(e_{i} + e_{ij})}{\sum_{k}^{}{e_{k}/(e_{i} + e_{ik})}}$
(D.3)

The total attractiveness of an ED is then estimated by Equation D.4,
which is then normalised as per Equation D.5.

$T_{j} = \ \sum_{i}^{}T_{ij}$ (D.4)

${\overline{T}}_{j} = \frac{{\overline{T}}_{j}}{\sum_{j}^{}{\overline{T}}_{j}}$
(D.5)

The implementation of the above formulations is done as following:

4.  We found the centroid of each ED as a reference point for population
    concentration in that ED.

5.  We assumed that the people in an ED would look for employment
    opportunities which are within 10 km from them, i.e.,
    $d_{ij} = 10km$.

6.  For each ED $i$ we created a 10 km radius around the centroid and
    located every other ED $j$ within that radius, which gave us the
    values of $e_{ij}$. We repeated this calculation for each ED in
    Jamaica to find $e_{ik}\forall k$, which then solved Equations D.3
    -- D.5.

7.  The process outlined in Step 3 allowed us to assign a weight to each
    ED in Jamaica to signify how likely the working population would be
    concentrated within that ED, and hence be likely to work in the
    buildings within that ED. We assumed that the spatial disaggregation
    of GDP of a sector would be proportional to the working population
    attracted to EDs, which meant that the GDP of a sector $s$ assumed
    to be concentrated within an ED $i$, ${GDP}_{si}$, was estimated as:

${GDP}_{si} = {\overline{T}}_{i}{GDP}_{s}$ (D.6)

8.  Within a given ED we then looked at all the buildings of that sector
    and assumed that GDP of sector $s$ assigned to building $b$ in ED
    $i$ was estimated from Equation D.7, where $a_{sib}$ denotes the
    area of the building tagged to sector $s$ and within ED $i$.

${GDP}_{sib} = \frac{a_{sib}}{\sum_{b}^{}a_{sib}}{GDP}_{si}$ (D.7)

## D.3 Disaggregating GDP to agriculture areas

1.  We estimate national scale GDP values from the steps outlined in
    Section D.1, for the agriculture sector and subsectors with codes A
    011-1 -- 011-8, A 12, A 14, A 20 shown in Table D-1.

2.  We use the FD and TNC land-use datasets to find all land-use types
    that can be mapped to the agriculture sector and subsector classes.
    For example, these datasets show the areas where plantation crops
    such as sugarcane and bananas are grown, which means we can map
    those areas to sector and subsector code A 011-1 and A 011-2 and
    later on disaggregate the GDP of these sectors to these areas.

3.  We use the IFPRI Map Spatial Production Allocation Model
    (MapSPAM)[^77] data estimates from 2010 to estimate the total value
    for different crops values produced in Jamaica. MapSPAM provides
    estimates production values in annual tonnages and US\$ for 42 crops
    globally at a 5km gridded resolution[^78]. Using this data, we map
    the crops to the agriculture subsectors for Jamaica and then
    estimate the total subsector production values per areas from 2010
    in US\$/m^2^ disaggregated at 5km gridded resolutions.

4.  We intersect the land-use datasets (step 2 above) with the 5km
    gridded MapSPAM data layers (step 3 above) and assign the production
    values per area of the specific sectors to the land use areas.

5.  We assume that the overall subsector GDP value (step 1 above) will
    be disaggregated to areas in proportion to the estimated production
    values in US\$/m^2^ estimated for these areas. Hence, we are using
    the 2010 estimate as weights to spatially disaggregate 2019 GDP
    estimates. We note that there is an assumption here that agriculture
    production patterns and intensity in 2019 will be similar to those
    in 2010, which is made in the absence of any other data on spatial
    agriculture production.

## D.4 Disaggregating GDP to mining and quarrying areas

1.  We estimate national scale GDP/day values from the steps outlined in
    Section D.1, for the mining and quarry subsectors with codes C 132
    and C 141.

2.  We use the FD and TNC land-use datasets to find all land-use types
    that can be mapped to the mining and quarrying area. For example,
    these datasets show the areas where mining and quarrying activities
    are taking place, which means we can identify all areas where the
    GDP associated with sectors C 132 and C 141 are coming from.

3.  We do not have any information on the mining and quarrying outputs
    (in tonnages or J\$) produced from specific locations of the mining
    and quarrying areas we have identified.

4.  For quarrying we assume that the quarrying outputs and hence GDP is
    disaggregated in proportion to the areas of the quarrying locations
    identified from the FD and TNC data. Since, we do not have any
    information on the intensity of output from each quarry in Jamaica.
    we are simply assuming the same output/area for each quarry.

5.  For mining we know the total tonnages of bauxite and alumina (C 132)
    being exported from ports in Jamaica as per statistics published by
    Port authority of Jamaica (PAJ)[^79]. We assume
    that Using the transport network of roads and railways we map the
    mines to the ports, which allows us to know which mines are
    contributing to the total export tonnages to specific ports.

6.  We assume that the contribution of a mine to a port is in proportion
    to its area, and hence we can estimate the total tonnage outputs of
    specific mines connected to each port in Jamaica.

7.  From the estimate of total tonnages allocate to each mine, we assume
    that the total GDP of the sector C 132 is disaggregated to the mines
    in the same proportion.

[^77]: <https://www.mapspam.info/>
[^78]: <https://www.mapspam.info/methodology/>
[^79]: <https://www.portjam.com/stat-report/Monthly_Statistical_Publication_February_2020.pdf>

## Contents

1. [Introduction](01-introduction.md)
2. [Methodology development and implementation steps](02-methodology.md)
3. [Model assumptions and data assembled for implementing J-SRAT](03-model-assumptions-data.md)

- [Appendix A: Vulnerability curves for infrastructure assets in Jamaica](appendix-a-vulnerability-curves.md)
- [Appendix B: Hazard models](appendix-b-hazard-models.md)
- [Appendix C: Infrastructure network flow models for failure analysis](appendix-c-network-flow-models.md)
- [Appendix D: Spatial disaggregation of economic activity at buildings and area levels](appendix-d-spatial-disaggregation.md)
