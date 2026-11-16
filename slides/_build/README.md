# Slide build sources

`build.py` turns every `*.src.html` here into a self-contained deck in `../dayN/`.

- `d3.min.js` – d3 v7.9 (ISC licence), inlined so decks run offline
- `gbsg2.json` – `tsize`, `age`, `pnodes` from the GBSG2 data (lifelines), used by the live simulations
- `logo_*.png` – downsized copies of the official GIMM logo files in `assets/gimm_logo/`

Edit the `.src.html`, then run `python slides/_build/build.py` from the repo root.
