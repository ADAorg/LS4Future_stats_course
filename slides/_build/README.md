# Slide build sources

Edit here, then run `python slides/_build/build.py` from the repo root. Every `DayN_*.src.html` becomes a
self-contained deck in `../dayN/` (one HTML file, works offline, safe to email or put on a USB stick).

| File | What it is | Inlined at |
|---|---|---|
| `DayN_Tk_topic.src.html` | One deck: its slides and slide-specific charts (`hooks`) | — |
| `shared/gimm-deck.css` | Shared GIMM deck styles | `/*CSS*/` |
| `shared/deck.js` | Shared engine: navigation, speaker notes, tooltips, simulation helpers, `startDeck(hooks)` | `/*DECKJS*/` |
| `shared/gbsg2.js` | GBSG2 variables and population summaries (Day 1 decks) | `/*GBSG2*/` |
| `gbsg2.json` | `tsize`, `age`, `pnodes` from GBSG2 (lifelines) | `__DATA__` |
| `lung.json` | NCCTG lung cancer data (lifelines): pat_karno, ecog, age, sex, wt_loss, time, status | `__JSON:lung__` |
| `airway_mv.json` | Airway RNA-seq (GSE52778), control samples: per-gene mean and variance after library-size scaling, genes with mean ≥ 1 | `__JSON:airway_mv__` |
| `airway_de.json` | PyDESeq2 0.6 results on airway, built in Python: paired (~celltype + dex) and unpaired (~dex) designs (log2FC, −log10 p, padj < 0.05), plus paired runs on each pair of cell lines for the winner's-curse slide | `__JSON:airway_de__` |
| `airway_top.json` | Top 10 paired genes (index → Ensembl ID, symbol where known) | `__JSON:airway_top__` |
| `cox3.json` | Cox models (lifelines) for GBSG2, lung and the HO1 simulated dataset (crude / + age / + all covariates), plus the simulated survival data itself | `__JSON:cox3__` |
| `d3.min.js` | d3 v7.9 (ISC licence) | `/*D3*/` |
| `logo_neg.png`, `logo_pos.png` | Downsized copies of the official GIMM logos | `__LOGO_NEG__`, `__LOGO_POS__` |

A change to anything in `shared/` reaches every deck on the next build.
Deck-specific styles can go in an extra `<style>` block in that deck's `.src.html`.
