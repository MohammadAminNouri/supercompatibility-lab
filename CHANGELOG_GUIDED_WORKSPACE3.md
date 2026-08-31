# Guided Workspace 3 revision — 2026-08-31b

## UX changes
- Workspace 3 now contains its own NiTi lattice inputs.
- One task is shown at a time: variants/operators, twins/shear-shear, or Interaction Work.
- Raw matrices and full 12x12/groupoid tables are advanced/audit views.
- Variant-pair explorer explains what an operator class means.
- Twin explorer is selection-first rather than table-first.
- IW accepts either B2 [h k l] loading or reconstructed parent B2 Bunge Euler angles.

## Scientific correction to Interaction Work
Previous integrated revision ranked B2->B19' formation distortions directly. This revision implements martensite reorientation from a declared initial state i:

F_i->j = F_j F_i^-1

epsilon_i->j = F_i->j - I

IW_i->j = sigma : epsilon_i->j

This follows the reorientation construction in Xiao, CT-reference-author & Loge, International Journal of Plasticity 159 (2022) 103468, and the IW definition reused in Scripta Materialia 256 (2025) 116433.

## Variant terminology
- V1..V12 in the operator task are correspondence variants from left cosets.
- D1..D12 in the IW task are explicit software-ordered distortion matrices.
- They are intentionally not silently equated by number. Exact D matrices are exportable.
