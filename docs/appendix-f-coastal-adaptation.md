# Appendix F: Coastal flood protection adaptation

Coastal flooding represents a significant hazard to Jamaica's infrastructure,
particularly given the concentration of critical assets in coastal zones and
projected sea-level rise under climate change scenarios. This appendix describes
the methodology developed for identifying and evaluating coastal flood
protection adaptation options, specifically engineered coastal defence
structures such as revetments, seawalls, and levees.

Unlike asset-specific adaptation measures (such as elevating individual
buildings or upgrading structural resilience), coastal flood protection provides
regional-scale defence by reducing flood extent and depth over contiguous
coastal areas. The methodology addresses the challenge of spatially delineating
protection zones, assigning infrastructure assets to these zones, apportioning
costs, and evaluating the collective benefit-cost ratio of protecting multiple
assets within each zone.

## Methodology

### Coastal flood protection zone delineation

The methodology identifies discrete coastal areas that would benefit from flood protection infrastructure and determines the corresponding coastline segments where defence structures would be constructed. This process uses the most severe coastal flood scenario (typically RCP 8.5, 100-year return period, year 2100) to define the maximum potential extent of coastal flooding that the adaptation measures should address.

**1. Initial flood area identification**

The DBSCAN (Density-Based Spatial Clustering of Applications with Noise) algorithm is applied to coastal flood depth raster maps to identify spatially distinct flooded areas. Areas exceeding a flood depth threshold (typically 0.1 m) are grouped into discrete clusters representing contiguous flood-affected zones. The algorithm parameters include:

- **Epsilon (eps)** -- Maximum distance between pixels to be considered part of the same cluster (typically 3 pixels)
- **Minimum points (minpts)** -- Minimum number of pixels required to form a cluster (typically 5)

For each identified cluster, a minimum enclosing polygon is generated representing the extent of flooding. This creates an initial set of discrete flood-affected areas along Jamaica's coastline.

**2. Coastline-based zone generation (SCAPE method)**

To create protection zones aligned with the coastline geometry, the methodology uses a coastline-based segmentation approach:

1. For each segment of the Jamaica coastline, a rectangular bounding box is constructed perpendicular to the coastline orientation, extending inland by a buffer distance (typically 2 km). This creates an initial set of coastal polygons aligned with the coastline.

2. The maximum flood depth within each coastal polygon is calculated from the flood hazard raster using zonal statistics.

3. Adjacent coastal polygons are grouped together based on flood height thresholds and coastline length constraints. Polygons are joined if:
   - Their flood heights are consistently above or below the threshold (typically 0.1 m)
   - The cumulative coastline length does not exceed the maximum segment length (typically 5 km)
   - Transitions between high and low flood heights trigger group finalization

4. Overlapping or adjacent groups are merged if their overlap exceeds a threshold (typically 50%), creating a consolidated set of coastal protection zones.

This process produces an initial set of flood protection zones that are aligned with the coastline geometry and respect length constraints for practical construction.

**3. Voronoi tessellation for length control**

Large flood protection zones are subdivided using Voronoi tessellation to ensure individual coastal defence segments remain within the maximum length constraint:

1. For each protection zone from step 2, multiple seed points are generated along the coastline in proportion to the zone's extent.

2. A Voronoi diagram is constructed from these seed points, creating Voronoi regions that subdivide the protection zones.

3. The flood polygons identified in step 1 are associated with each Voronoi region by spatial intersection, determining which flooding each Voronoi region would defend against.

4. For each Voronoi region, the coastline segment within that region is extracted, creating individual coastal defence segments that collectively protect the flood-affected areas.

This tessellation ensures that even large continuous flood zones are divided into manageable coastal segments for defence construction, with each segment defending a portion of the overall flood area.

**4. Protection area creation and buffering**

For each Voronoi-derived coastline segment, a flood protection area polygon is generated:

1. The coastline segment is buffered inland, starting with a default buffer distance (typically 400 m).

2. The buffer is iteratively expanded until it fully encloses the flood polygon associated with that Voronoi region.

3. The buffered area is clipped to the Voronoi region boundary to prevent overlap between adjacent protection areas.

This creates discrete protection area polygons, each associated with a specific coastline segment, that fully encompass the flooding they are designed to defend against.

**5. Refinement and consolidation**

The protection areas and coastline segments undergo final refinement:

1. Multi-part polygons (resulting from irregular flood geometries) are split into separate components, with small isolated sections identified for potential removal or merger.

2. Adjacent protection areas along the coastline are evaluated for merging if:
   - Their combined coastline length remains under the maximum limit
   - They share significant boundary contact
   - Merging would create more efficient protection schemes

3. The maximum flood depth within each final protection area is calculated from the flood hazard raster and assigned as the design height for the associated coastal defence segment.

4. Coastline segments are extracted by intersecting the full Jamaica coastline with each final protection area, creating linestring geometries that define where defence structures would be built.

**Output structure**

The process produces two primary outputs:

- **Coastal protection feature areas** -- Polygons representing inland zones that would be defended by coastal structures, each with an associated maximum flood height
- **Coastal defence segments** -- Linestrings along the shoreline where structures would be built, with associated design heights and lengths

These outputs form the spatial basis for subsequent asset mapping and benefit-cost analysis.

### Asset-to-protection zone mapping

Infrastructure assets (roads, electricity substations, water treatment facilities, etc.) within coastal areas must be mapped to the flood protection zones that would defend them. This mapping enables:

1. Identification of which assets benefit from each protection option
2. Calculation of avoided damages and losses for each protection zone
3. Apportionment of protection costs among benefiting assets

The mapping process involves:

**1. Spatial intersection**

For each infrastructure asset (point, line, or polygon), a spatial intersection is performed with all coastal flood protection areas across multiple climate scenarios (different RCPs, return periods, and time epochs). An asset is considered "protected" by a given coastal defence if any part of it falls within that defence's protection area.

**2. Protection height assignment**

For each asset-protection zone pairing, the flood height that the coastal defence structure would need to achieve is recorded. This is derived from the maximum flood depth within the protection area for the given climate scenario.

**3. Multi-scenario tracking**

Assets may be mapped to different protection zones under different climate scenarios (e.g., a modest coastal protection might suffice under RCP 4.5 but a more extensive zone might be needed under RCP 8.5). The mapping records these relationships across all scenarios.

### Cost calculation and apportionment

The cost of implementing a coastal flood protection scheme includes construction and ongoing maintenance. For Jamaica, cost estimates are based on:

- **Construction costs** -- Unit costs per linear meter of coastal defence structure (seawall, revetment, etc.) as a function of height
- **Maintenance costs** -- Annual routine maintenance and periodic major rehabilitation

**Cost apportionment among assets**

When multiple infrastructure assets benefit from a single coastal protection zone, the total protection cost must be apportioned. The approach used is:

1. Calculate the total rehabilitation cost of all infrastructure assets within the protection zone
2. For each individual asset, assign a fraction of the total coastal protection cost proportional to that asset's rehabilitation cost relative to the total

Mathematically, for an asset $i$ within protection zone $z$:

$$\text{Apportioned Cost}_i = \text{Total Protection Cost}_z \times \frac{\text{Asset Rehab Cost}_i}{\sum_{j \in z} \text{Asset Rehab Cost}_j}$$

This ensures that higher-value assets bear a proportionally larger share of the protection cost, reflecting the principle that protection investment should be commensurate with the value at risk.

### Benefit-cost analysis

The effectiveness of each coastal protection option is evaluated through a cost-benefit analysis, following the general adaptation assessment framework described in Section 2.2.

**Benefits calculation**

For each asset protected by a coastal defence option:

1. **Without adaptation** -- Calculate EAD and EAEL under coastal flood hazard for all return periods, climate scenarios, and time epochs
2. **With adaptation** -- Recalculate EAD and EAEL assuming the coastal protection eliminates flooding up to the design flood height. Assets in areas with flood depths below the protection height experience zero coastal flood damage; those in areas exceeding the protection height still incur damages from the residual flooding.
3. **Avoided risks** -- The difference between scenarios (1) and (2) represents the avoided EAD and EAEL, which constitute the benefits of the adaptation option

**Aggregation levels**

Benefit-cost ratios are calculated at two levels:

1. **Per-asset BCR** -- For each individual infrastructure asset, comparing its apportioned share of the coastal protection cost against the avoided damages and losses it (doesn't) experience. This allows identification of which specific assets benefit most from coastal protection.

2. **Per-protection-zone BCR** -- For each coastal protection feature area, aggregating the benefits (avoided damages and losses) across all protected assets and comparing against the total protection cost. This provides an overall assessment of whether the protection scheme is economically justified for that coastal segment.

**Multi-hazard consideration**

Assets may be exposed to multiple hazard types beyond coastal flooding (fluvial flooding, cyclones, etc.). The coastal adaptation benefit-cost analysis focuses exclusively on avoided coastal flood risks. Other hazards and their adaptation options are assessed independently, though in practice, decision-makers may consider the combined benefits of assets that warrant protection under multiple hazard scenarios.

## Application and outputs

The coastal flood protection analysis produces several key outputs:

**Spatial datasets**

- **Coastal protection feature areas** -- Polygons representing inland zones defended by coastal structures
- **Coastal defence segments** -- Linestrings along the shoreline where structures would be built, with associated design heights
- **Asset-protection mappings** -- Tables linking infrastructure assets to their protecting coastal defence zones across climate scenarios

**Economic analysis results**

- **Per-asset adaptation summaries** -- For each infrastructure asset within a coastal zone: apportioned protection cost, avoided EAD, avoided EAEL, NPV costs, NPV benefits, and BCR
- **Per-protection-zone summaries** -- For each coastal defence segment: total construction cost, aggregated avoided damages across all protected assets, total NPV benefits, and BCR
- **Return period damages and losses** -- For each protection zone, the avoided damages and losses at each return period under protected versus unprotected scenarios

**Prioritization and planning**

The BCR values enable prioritization of coastal protection investments:

- Protection zones with BCR ≥ 1 indicate economically justified interventions
- Ranking zones by BCR identifies highest-priority locations
- Assets with high avoided losses but low BCRs may indicate cases where the protection cost exceeds the benefits, suggesting alternative adaptation strategies might be more appropriate

The spatial outputs also support planning by:

- Identifying coastal segments where protection is most needed
- Informing required design heights for coastal structures
- Highlighting co-benefits where a single coastal defence protects multiple infrastructure sectors (energy, transport, water)

## Limitations and considerations

Several important considerations apply to this methodology:

**Simplicity of flood protection model**

We assume that building a coastal defence segment between the sea and an area deemed at flood risk will eliminate that risk. Without an inundation model and an accurate map of terrain, this is a large leap. For example, unless all neighbouring zones are built, flooding may come round a single wall section and inundate the area behind. However, we think this simple method is justified for a 'first look' at potential protection zones.

**Design standard trade-offs**

Coastal defences designed to the "worst case" scenario (e.g., RCP 8.5, 100-year flood, year 2100) may be overdesigned for moderate scenarios but provide robust protection under deep uncertainty about future climate. Alternative approaches might consider incremental adaptation pathways with lower initial costs.

**Residual risk**

Coastal protection reduces but does not eliminate flood risk. Assets in protected zones may still experience flooding from:

- Events exceeding the design standard
- Defence structure failure or breach
- Fluvial or pluvial flooding originating from inland areas

**Nature-based solutions**

This methodology focuses on engineered coastal defences (hard infrastructure). Nature-based solutions such as mangrove restoration, coral reef protection, or beach nourishment may offer alternative or complementary approaches with different cost profiles and co-benefits (ecosystem services, tourism value). These are considered separately in the J-SRAT framework.

**Cost uncertainty**

Construction cost estimates are subject to significant uncertainty, particularly for structures in remote or difficult terrain. Sensitivity analysis across cost ranges is recommended when using BCR results for decision-making.
