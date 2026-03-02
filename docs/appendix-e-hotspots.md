# Appendix E: Spatial hotspots analysis

While asset-level risk metrics identify individual critical 
infrastructure elements, spatial hotspots analysis provides a 
complementary regional perspective that identifies geographic areas where 
infrastructure risks and vulnerabilities are concentrated. This analysis 
produces spatially aggregated risk maps that support strategic planning 
and investment prioritisation at regional and national scales.

The hotspots analysis framework combines three types of spatial metrics:

1. **Exposure value hotspots** -- Geographic distribution of 
   infrastructure asset replacement costs, aggregated by grid cell
2. **Risk hotspots** -- Spatial distribution of Expected Annual Damages 
   (EAD) across hazard types and sectors
3. **Economic loss hotspots** -- Regional criticality based on the 
   wider economic consequences of infrastructure disruption within each 
   area

## E.1 Methodology

The methodology for generating these hotspots comprises the following steps:

### Grid definition

A regular raster grid is defined to cover Jamaica's extent, with a configurable 
cell resolution (typically 1 km × 1 km). This grid serves as the spatial unit 
for aggregating infrastructure metrics.

### Asset spatial disaggregation

Infrastructure assets (point, line and polygon geometries) are spatially 
intersected with the hotspots grid. Line and polygon assets that span multiple 
grid cells are split at cell boundaries, with geometric properties (length, 
area) recalculated for each split segment. Each split retains the original 
asset's attributes, including rehabilitation costs per unit dimension.

### Exposure value calculation

For each grid cell, the total exposure value is calculated by summing the 
rehabilitation costs of all asset segments within that cell. For point assets 
(e.g., substations, treatment plants), the full asset cost is assigned to the 
containing cell. For linear assets (e.g., roads, transmission lines), the cost 
is calculated as the product of the split segment length and the unit cost 
(J\$/m). For polygon assets (e.g., ports, airports), the cost is the product of 
the split area and unit cost (J\$/m²). These cell-level exposure values are 
aggregated by sector (energy, transport, water) and across all sectors to 
produce exposure hotspot maps.

### Risk hotspots from Expected Annual Damages

Building on the direct damage assessment methodology described in Section 2.2, 
EAD values are calculated for each asset split segment across all hazard types, 
return periods, and climate scenarios. These asset-level EAD estimates are then 
spatially aggregated to the hotspots grid by summing the EAD values of all asset 
splits within each cell. Risk hotspots are produced for individual hazard types:

- coastal flooding
- fluvial flooding
- surface water flooding
- combined flooding (sum of all flood types)
- tropical cyclones

Along with total multi-hazard risk. These metrics are calculated per sector and 
for all infrastructure combined.

### Transport economic loss hotspots

For the transport sector, a network disruption analysis quantifies the wider 
economic consequences of losing all road and rail infrastructure within each 
grid cell. This process involves:

1. For each grid cell, identifying all road and rail edges that intersect with 
   that cell
2. Removing these edges from the multi-modal transport network (as described in 
   Appendix C) and re-optimising all origin-destination flows
3. Calculating the economic impacts of this disruption on both trade flows and 
   labour flows:

   - **Trade losses** -- Measured as the increased cost of rerouting freight 
     flows (or complete loss of access to isolated origins or destinations), 
     quantified in J\$/day
   - **Labour losses** -- Measured as the cost of additional travel time for 
     commuters (calculated using the hourly labour cost and increased journey 
     time), or the GDP loss from workers unable to reach their place of 
     employment, quantified in J\$/day

The total economic loss for each cell is the sum of trade rerouting losses, 
trade isolation losses, labour time losses, and labour GDP losses. This analysis 
is performed independently for each cell using the nominal (baseline) transport 
network and flow patterns.

### Spatial smoothing

To enhance the visual interpretability of hotspot maps and account for the 
spatial continuity of infrastructure risks, an optional Gaussian kernel density 
estimation (KDE) smoothing is applied to the gridded results. This 
quantity-preserving smoothing process redistributes the total value in each 
coarse grid cell to a finer resolution output grid (typically 250 m × 250 m) 
using a Gaussian kernel with configurable bandwidth (typically 1 km). The KDE 
approach ensures that the total summed value across Jamaica remains constant 
while creating smoother, more continuous spatial patterns that better represent 
the regional nature of infrastructure risk.

## E.2 Application

The hotspots analysis provides decision-makers with spatially explicit 
information about where infrastructure risks are concentrated, complementing the 
asset-level criticality assessments. These maps can guide strategic investment 
in resilience measures, inform land-use planning to avoid high-risk areas, and 
support the prioritisation of adaptation interventions at regional scales.
