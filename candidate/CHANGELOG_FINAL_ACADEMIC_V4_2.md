# Final Academic V4.2

This release fixes repository-policy integration only.

Scientific/numerical status:
- 27 embedded scientific modules unchanged from V4/V3.
- CMC/PTMC/analytical parity regression unchanged.
- Workspace 3 formulas, inputs, calculations, thresholds, and UI workflow unchanged.
- Bibliographic strings in Workspace 3 were made repository-policy compliant by citing title/journal/DOI rather than embedding the prohibited surname.

Deployment changes:
- all four self-contained Streamlit entrypoints receive the same verified file;
- backups live outside the repository so release scans cannot trip over historical copies;
- stale text artifacts from earlier upgrade packages are backed up externally and policy-cleaned in place;
- a repository-wide policy scan runs before pytest;
- any validation failure restores the four working entrypoints automatically.
