# Parent ↔ daughter reconstruction module

This document describes the grain-level orientation reconstruction workbench and the conventions needed to use it reproducibly.

## Scope

The module solves two related problems:

1. **Forward variant generation:** given a parent orientation and a crystallographic orientation relationship (OR), enumerate the symmetry-distinct daughter orientations.
2. **Inverse parent reconstruction:** given measured daughter grain orientations plus a neighborhood graph, infer groups of daughter grains likely inherited from the same parent and estimate the parent orientation.

It is intentionally independent of EBSD vendor file formats. Import/segmentation of raw EBSD pixels and grain construction should be performed before this module. The input is grain-level orientation data.

## Internal orientation convention

The internal orientation matrix \(G\in SO(3)\) maps a vector written in **crystal coordinates** into the **specimen frame**.

Uploaded data can be declared as either:

- crystal → specimen, used directly; or
- specimen → crystal, explicitly transposed on import.

This choice is not cosmetic. Reversing an orientation matrix changes OR composition and can invalidate a reconstruction even if the Euler angles look plausible. The source software/export convention should therefore be recorded in any manuscript or archived analysis.

Bunge Euler inputs are in degrees and quaternions use scalar-first columns `qw,qx,qy,qz`.

## OR representation

The OR is stored as a proper child→parent crystal-frame rotation \(R_{cp}\). An ideal daughter variant generated from parent orientation \(G_p\) is

\[
G_c=G_p S_p R_{cp},
\]

where \(S_p\) is a proper parent symmetry. Equivalent daughter orientations are removed modulo the daughter point group.

For one measured daughter orientation \(G_c\), the inverse problem generates all parent candidates consistent with the OR and daughter symmetry. Candidate parents are deduplicated modulo parent symmetry.

## Built-in OR families

| OR family | Parent → daughter | Plane relation | Direction relation | Tested ideal variant count |
|---|---|---|---|---:|
| Kurdjumov–Sachs | FCC → BCC | `{111}p || {110}d` | `<110>p || <111>d` | 24 |
| Nishiyama–Wassermann | FCC → BCC | `{111}p || {110}d` | `<011>p || <001>d` | 12 |
| Bain | FCC → BCC | `{100}p || {100}d` | `<100>p || <110>d` | 3 |
| Pitsch | FCC → BCC | `{100}p || {110}d` | `<110>p || <111>d` | 12 |
| Burgers | BCC → HCP | `{110}p || (0001)d` | `<111>p || <11-20>d` | 12 |

The code tests the rotation matrices against the plane/direction parallelisms as well as the symmetry-generated variant counts.

## Neighborhood data

The preferred adjacency input is a table of actual neighboring daughter grains from a segmented EBSD map:

```text
grain_id_1,grain_id_2
17,18
18,25
...
```

A centroid k-nearest-neighbor graph can be generated when `x,y` columns are present, but it is explicitly marked **approximate** because nearest centroids need not share a physical grain boundary.

## Reconstruction families

### Neighbor voting

For grain \(i\), candidate parent \(P_{ik}\) receives support from each adjacent grain \(j\) according to the best parent-symmetry-reduced angular fit to one of \(j\)'s candidates. The highest-support candidate is selected locally; spatially adjacent selected parents are merged within the user-specified parent merge tolerance.

Use this method as a simple and interpretable baseline.

### Grain graph + MCL

A graph node represents a daughter grain. Edge weights decrease with the best common-parent OR misfit. The graph is partitioned by Markov clustering. The exposed inflation parameter controls cluster granularity and must be reported if the result is used in a publication.

This reflects the established grain-graph reconstruction strategy in which OR-consistent neighboring daughter grains are clustered into prior parents.

### Candidate-level variant graph

Each daughter grain contributes one node per possible parent candidate. Compatibility between candidate nodes on neighboring grains is converted into edge support. Probability-like support is propagated iteratively before the best parent candidate is selected.

The important conceptual difference is that second- and third-best candidate hypotheses are retained during graph propagation rather than being discarded before clustering.

### Nucleation + growth

A strict angular threshold is used to identify locally well-supported parent seeds. Unassigned neighboring grains are then incorporated with a looser growth tolerance. This mimics established reconstruction workflows that separate confident nucleation from more permissive growth.

### Operator / groupoid consistency

The OR and parent point group generate a theoretical set of daughter–daughter operators. A measured neighboring pair first has to match one of these operators modulo daughter symmetry. Candidate parent orientations are then required to agree before the edge can connect a reconstructed parent region.

This method is useful as a crystallographically explicit counterpoint to purely weighted graph clustering.

## OR refinement

The OR refinement routine accepts a physically meaningful starting OR and optimizes a bounded 3-component rotation-vector correction. For every trial OR it reconstructs candidate-parent sets and evaluates neighboring grains by their best shared-parent mismatch. A pseudo-Huber loss reduces the influence of outlier boundaries. The numerical search is a deterministic coarse-to-fine coordinate pattern search, chosen because candidate switching makes the objective non-smooth and because deterministic bounded evaluations are preferable for reproducible CI and interactive use.

Outputs include:

- initial objective;
- final objective;
- correction angle in degrees;
- refined OR matrix.

A refinement that improves the objective does **not** by itself prove that the recovered OR is physically correct. OR refinement should be validated against retained-parent measurements, known pole parallelisms, independent software or a sufficiently rich set of daughter variants.

## Output fields

The grain table contains:

- `reconstructed_parent_id` — cluster label;
- `variant_candidate_id` — selected parent-candidate / OR variant index under the internal enumeration;
- `fit_deg` — parent-symmetry-reduced angular fit to the reconstructed parent;
- `confidence` — heuristic 0–1 support score based on absolute fit and separation from the next candidate; **not a calibrated probability**;
- parent Bunge Euler angles;
- parent scalar-first quaternion.

Method-comparison exports also report runtime, parent count, mean/median fit and mean support score. Synthetic validation runs add permutation-invariant clustering accuracy against known parent labels.

## Recommended reporting

At minimum, report:

- input OR family and any refined correction;
- parent and daughter point groups;
- orientation matrix direction convention;
- grain segmentation source and adjacency construction;
- method and all angular/graph parameters;
- number of reconstructed parents;
- fit distribution, not only the mean;
- support-score distribution, clearly labelled heuristic;
- disagreement between at least two reconstruction methods near important boundaries;
- any retained-parent or independent validation used.

## Scale and performance

The implementation prioritizes transparent equations and testability over specialized sparse-graph optimization. Candidate-level graph methods become expensive as the number of grains and variants grows. Large production EBSD maps containing many thousands of grains are better processed with an optimized reconstruction platform, then imported here for cross-checking or region-level analysis.
