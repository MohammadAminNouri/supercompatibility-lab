# Martensite Research Suite

A clean-room, paper-oriented crystallography application intended to go beyond a simple variant generator.

## What is implemented now

- Generic parent/daughter proper rotational symmetries
- Literature-grounded NiTi B2→B19′ natural AQ OR preset
- Common-subgroup calculation
- Orientational variants from cosets
- Double-coset operator classes
- Variant-to-variant operator matrix
- Multivalued operator composition table
- Pole-figure-style direction plotting
- Interaction-work calculation `IW = sigma:epsilon`
- Uniaxial loading sweeps and IW-max direction search
- Editable literature NiTi transformation-gradient preset
- Right-stretch eigenvalues and `|lambda2-1|` diagnostic
- Rank-one distance diagnostic
- Known-truth reconstruction validation metrics
- Paper-ready CSV/JSON/equation/reference evidence ZIP
- Explicit model scope and provenance

## Important scope limits

This first public build does **not** claim to implement the complete 2026 Correspondence-Theory CMC/SMC/shear-shear supercompatibility derivation. It deliberately labels those as future/advanced modules rather than fabricating equations.

Likewise, it does not claim that an IW maximum alone predicts kinetics. The IW module calculates the mechanical-work term and ranking under explicitly declared stress/strain inputs.

## Run

```bash
python -m pip install -r requirements.txt
streamlit run app.py
```

## Test

```bash
pytest -q
```
