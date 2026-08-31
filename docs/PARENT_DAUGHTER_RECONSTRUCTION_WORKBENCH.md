# Parent ↔ daughter reconstruction workbench

## Scientific input contract

The reconstruction engines operate on **daughter grains**, not on an unexplained collection of Euler angles.  A defensible reconstruction requires:

1. one measured daughter orientation per segmented daughter grain;
2. the parent and daughter proper-rotation symmetry groups;
3. a parent/daughter orientation relationship (OR), supplied from literature, explicit plane/direction parallelisms, or an explicit rotation matrix;
4. a graph describing which daughter grains physically neighbor each other;
5. explicit angular/graph thresholds for the selected reconstruction method.

`grain_id` is a unique daughter-grain identifier.  Orientations may be supplied as Bunge Euler angles `phi1_deg, Phi_deg, phi2_deg` in degrees, or as unit quaternions `qw,qx,qy,qz`.  The software's internal orientation convention is crystal → specimen and an alternative specimen → crystal input is transposed explicitly.

## Input routes

The Streamlit workbench offers four main routes:

- built-in two-parent validation data;
- a manual grain table and manual/measured/approximate adjacency;
- a pre-segmented grain CSV/TSV/TXT;
- raw `.ang` or `.ctf` EBSD points followed by an explicit in-app spatial/misorientation segmentation step.

The raw-point segmenter is intended for modest maps, method education and auditable first-pass studies.  It is not presented as a replacement for mature EBSD segmentation in MTEX, OIM Analysis, AZtecCrystal or comparable tools.  Production-scale studies should generally import a pre-segmented grain table and measured grain-boundary adjacency.

## Orientation relationship

Three OR input modes are available:

- named literature preset;
- one parent plane/direction parallel to one daughter plane/direction;
- explicit 3×3 daughter→parent OR matrix `R_cp`.

Optional bounded OR refinement adjusts a supplied physical OR only within a user-reported angular bound.  It does not discover an unconstrained OR from scratch.

## Reconstruction families

The same dataset can be processed by one or several methods simultaneously:

- **Neighbor voting** — local candidate-parent support from adjacent grains.
- **Grain graph + Markov clustering** — one graph node per daughter grain, with OR-consistency edge weights.
- **Variant graph probability propagation** — candidate-parent nodes are retained explicitly before selecting parent assignments.
- **Nucleation + growth** — strict parent seeds followed by looser spatial growth.
- **Operator/groupoid consistency** — daughter/daughter transformation operators plus common-parent consistency.

Cross-method partition agreement is reported with the Adjusted Rand Index (ARI).  ARI compares clustering only; it does not prove that the reconstructed parent orientation is physically correct.

## Output hierarchy

### Parent summary

One row per reconstructed parent contains:

- reconstructed parent ID;
- number of supporting daughter grains;
- supporting area and area fraction when area is supplied;
- mean, median, 95th-percentile and maximum OR angular residual;
- mean/minimum heuristic support score;
- number of weak/ambiguous daughter assignments;
- parent Bunge Euler orientation;
- parent unit quaternion;
- a human-readable heuristic quality flag.

### Daughter assignment table

One row per daughter grain contains:

- daughter `grain_id`;
- reconstructed `parent_id`;
- selected candidate-variant ID;
- best OR angular residual in degrees;
- heuristic 0–1 support score;
- quality flag;
- available map coordinates, area and phase metadata;
- reconstructed parent orientation.

The support score is **not a calibrated probability**.  It combines absolute angular fit and separation from alternative candidates in the underlying numerical engine.

### Diagnostic plots

The workbench provides a reconstructed parent map when x/y are available, OR-fit distributions, support-score distributions and cross-method comparison tables.

## Forward problem

Known parent orientations can be entered manually or uploaded in bulk.  The app generates every symmetry-distinct daughter orientation allowed by the selected OR and exports both Euler angles and quaternions.

This is a **crystallographic variant prediction** only.  It does not determine nucleation site, selected variant, volume fraction, morphology or spatial arrangement without additional thermomechanical/microstructural physics.

## Minimum reportable Methods information

A publication using the reconstruction workbench should report at least:

- parent and daughter phases and symmetries;
- source/definition of the OR and any OR refinement;
- orientation representation and convention;
- daughter-grain segmentation method and threshold;
- adjacency source (measured boundary topology vs approximate k-NN);
- reconstruction method(s);
- all method-specific angular and graph parameters;
- number of daughter grains and adjacency edges;
- reconstructed parent count;
- OR residual statistics and support-score definition;
- retained-parent / known-truth validation if available;
- software build/commit and exported reconstruction tables.
