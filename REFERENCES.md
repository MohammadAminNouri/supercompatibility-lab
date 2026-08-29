# Primary literature map

This file records the literature used to define equations, benchmarks and research-extension modules. It is deliberately organized by **what the source supports** so that software features can be traced to evidence.

## Core compatibility theory

### Full cofactor conditions

**Study of the cofactor conditions: Conditions of supercompatibility between phases**  
DOI: <https://doi.org/10.1016/j.jmps.2013.08.004>

Software use:

- CC1 \(\lambda_2=1\)
- full CC2 cofactor equation
- CC3 inequality
- conventional Type-I/Type-II rank-one domain formulas
- simplified Type-I/Type-II CC2 checks under CC1

### Metric/correspondence compatibility formulation

**Compatibilities and supercompatibility conditions in shape memory alloys determined from correspondence, metrics and symmetries**  
DOI: <https://doi.org/10.1016/j.actamat.2026.122399>

Software use:

- built-in B2→B19′ correspondence convention
- CMC matrix and degeneracy interpretation
- SMC matrix and A/M shear extraction
- metric transformation-twin equations
- shear/shear mismatch \(\varepsilon\)
- B2/B19′ numerical benchmark and rounded target geometries
- distinction between demonstrated implication and unproved full converse equivalence

## Compatibility, hysteresis and reversibility benchmarks

### λ2-guided alloy tuning

**Energy barriers and hysteresis in martensitic phase transformations**  
DOI: <https://doi.org/10.1016/j.actamat.2009.05.034>

Database use: composition-tuned TiNiX benchmarks and strong hysteresis reduction near \(\lambda_2=1\).

### Near-cofactor cycling/reversibility benchmark

**Enhanced reversibility and unusual microstructure of a phase-transforming material**  
DOI: <https://doi.org/10.1038/nature12532>

Database use: Zn–Au–Cu high-reversibility benchmark, large transformation strain and highly stable transformation temperature during repeated cycling.

### Phase-engineering review

**Phase engineering and supercompatibility of shape memory alloys**  
DOI: <https://doi.org/10.1016/j.mattod.2017.10.002>

Software use: motivates keeping crystallographic compatibility, grain size and coherent precipitates as separate but jointly relevant research descriptors.

## Microstructure and fatigue

### Ultralow-fatigue precipitate/grain benchmark

**Ultralow-fatigue shape memory alloy films**  
DOI: <https://doi.org/10.1126/science.1261164>

Database use: Ti–Ni–Cu film benchmark with very high transformation-cycle durability and Ti2Cu precipitate/nanostructure context.

### Functional-fatigue mechanisms

**Functional fatigue and restoration in superelastic NiTi shape-memory alloys**  
DOI: <https://doi.org/10.1016/j.ijplas.2025.104483>

Database use: dislocations and retained martensite are tracked as mechanism-relevant metadata; no universal fatigue predictor is inferred.

### Grain-size distribution

**Role of grain size distribution in regulating the phase transformation behavior of polycrystalline NiTiCu shape memory alloys**  
DOI: <https://doi.org/10.1016/j.mtnano.2026.100893>

Database use: motivates recording distribution shape rather than only a mean grain size; the cited work is simulation-based and is labeled accordingly.

### Grain-size engineered elastocaloric cycling

**Achieving large room-temperature elastocaloric effect and ultrahigh cyclic stability by grain size engineering**  
DOI: <https://doi.org/10.1016/j.actamat.2026.122332>

Database use: recent experimental grain-size / \(\Delta T_{ad}\) / cyclic-stability benchmark.

## Elastocaloric and device performance

### Bulk Ni–Ti–Cu–Co

**Large room-temperature elastocaloric effect in a bulk polycrystalline Ni-Ti-Cu-Co alloy with low isothermal stress hysteresis**  
DOI: <https://doi.org/10.1016/j.apmt.2020.100844>

Database use: \(\Delta T_{ad}\), stress hysteresis, recoverable strain and COP context.

### Wide-temperature Ni–Ti–Cu–Fe

**Small stress-hysteresis in a nanocrystalline TiNiCuFe alloy for elastocaloric applications over wide temperature window**  
DOI: <https://doi.org/10.1016/j.jallcom.2022.167195>

Database use: wide superelastic temperature window, stress hysteresis, elastocaloric estimate and compatibility context.

### Two-step B2→R→B19′ transformation

**Low-fatigue large elastocaloric effect in NiTi shape memory alloy enabled by two-step transition**  
DOI: <https://doi.org/10.1016/j.scriptamat.2024.116239>

Software use: motivates the generic multi-step transformation-chain workbench.

### Elastocaloric device benchmark

**A zero-degradation elastocaloric cooling device using fatigue-resistant refrigerant**  
DOI: <https://doi.org/10.1016/j.joule.2026.102627>

Database use: device-level cooling power, cycle count, material temperature change and system temperature-span context.

## Data-driven design

### Composition→lattice ML precedent

**Predictions of Lattice Parameters in NiTi High-Entropy Shape-Memory Alloys Using Different Machine Learning Models**  
DOI: <https://doi.org/10.3390/ma17194754>

Software use: motivates user-trained composition/processing→lattice regression rather than hard-coding a fabricated universal model.

### Bayesian alloy/process optimization

**Design of high-temperature NiTiCuHf shape memory alloys with minimum thermal hysteresis using Bayesian optimization**  
DOI: <https://doi.org/10.1016/j.actamat.2024.120651>

Software use: precedent for inverse design over composition and processing variables; the current app limits itself to transparent lattice-space optimization unless the user supplies training data.

### Data-driven hysteresis prediction

**Data driven prediction of thermal hysteresis in NiTi-based high entropy shape memory alloys**  
DOI: <https://doi.org/10.1016/j.matlet.2025.140002>

Database use: current trend toward larger SMA datasets and interpretable model validation.

## Research frontier

### Extreme compatibility in compound domains

**Revisiting the cofactor conditions: Elimination of transition layers in compound domains**  
DOI: <https://doi.org/10.1016/j.jmps.2025.106409>

Software use:

- stores exact literature reference tensors for the stated cubic→monoclinic-II class;
- motivates stretch-variant commutator diagnostics;
- most importantly, motivates an **applicability guard** so transformation-class-specific theorems are not used outside their assumptions.

## Curated table

Machine-readable records are stored in [`../data/literature.csv`](../data/literature.csv). Empty cells mean that the curated evidence used for that row did not provide the quantity, not that the quantity is zero.

## Parent/daughter orientation reconstruction

### Groupoid / operator reconstruction foundations
**Groupoid of orientational variants** (Acta Crystallographica A, 2006)  
DOI: `10.1107/S010876730503686X`  
Supports the use of variants as cosets, inter-variant operators as double-coset classes, groupoid composition, and determination of possible parent crystals from daughter variants.

### Automated operator-based reconstruction
**ARPGE: a computer program to automatically reconstruct the parent grains from electron backscatter diffraction data** (Journal of Applied Crystallography, 2007)  
DOI: `10.1107/S0021889807048777`  
Supports operator-based neighboring-daughter checks, nucleation/growth reconstruction, parent-orientation recovery and variant/operator statistics.

### Least-squares OR-based parent reconstruction
**Mapping the parent austenite orientation reconstructed from the orientation of martensite by EBSD and its application to ausformed martensite** (Acta Materialia, 2010)  
DOI: `10.1016/j.actamat.2010.08.001`  
Supports least-squares fitting of parent orientation from measured daughter orientations under a specified/refined OR.

### Grain-graph and modern reconstruction workflow
**Parent grain reconstruction from partially or fully transformed microstructures in MTEX** (Journal of Applied Crystallography, 2022)  
DOI: `10.1107/S1600576721011560`  
Supports grain-graph reconstruction, OR refinement, parent voting and modern grain-level reconstruction workflows.

### Variant graph reconstruction
**The variant graph approach to improved parent grain reconstruction** (2022)  
See the MTEX `parentGrainReconstructor` documentation and associated paper/preprint.  
Supports retaining candidate parent variants as graph states and propagating candidate support rather than collapsing each child grain to one parent candidate immediately.

### Standard orientation relationships
The built-in KS, NW, Bain, Pitsch and Burgers presets use the conventional plane/direction parallelisms widely used in steels and titanium alloys. Their symmetry-generated ideal variant counts are software invariants and are tested explicitly.

## Compatibility methods beyond CC1–CC3

### Hadamard / general interface compatibility
**Compatibility conditions for microstructures and the austenite–martensite transition** (Materials Science and Engineering A, 1999)  
DOI: `10.1016/S0921-5093(99)00377-9`  
Supports use of the Hadamard jump condition for coherent planar interfaces and more general martensitic microstructures.

### Classical cofactor conditions
**Study of the cofactor conditions: Conditions of supercompatibility between phases** (Journal of the Mechanics and Physics of Solids, 2013)  
DOI: `10.1016/j.jmps.2013.08.004`  
Supports CC1–CC3, Type-I/Type-II specialization, and the all-twin-volume-fraction compatibility implication.

### Further supercompatibility and physically motivated distances
**On the cofactor conditions and further conditions of supercompatibility between phases** (Journal of the Mechanics and Physics of Solids, 2019)  
DOI: `10.1016/j.jmps.2018.08.012`  
Supports physically motivated distance concepts and additional martensitic-laminate supercompatibility structures. The current app cites this work but does not claim to implement every distance or star-twin theorem from it.

### Triplet condition
**Triplet condition: A new condition of supercompatibility between martensitic phases** (Journal of the Mechanics and Physics of Solids, 2022)  
DOI: `10.1016/j.jmps.2022.105050`  
Supports the triplet condition as a distinct three-variant supercompatibility principle.

### Triplet corrigendum
**Corrigendum to “Triplet condition: A new condition of supercompatibility between martensitic phases”** (Journal of the Mechanics and Physics of Solids, 2023)  
DOI: `10.1016/j.jmps.2023.105277`  
Must be consulted together with the primary triplet-condition paper because it corrects typographical errors in stated theoretical results.

### Compound-domain extreme compatibility
**Revisiting the cofactor conditions: Elimination of transition layers in compound domains** (Journal of the Mechanics and Physics of Solids, 2026)  
DOI: `10.1016/j.jmps.2025.106409`  
Supports transformation-class-specific extreme compatibility for compound domains. The application guard in this repository prevents these exact target structures from being used outside their stated class.
