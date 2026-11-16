# Practical Statistics for Biomedical Researchers

A 3-day introductory workshop on statistics and data analysis for molecular and clinical researchers, run in Google Colab with AI assistants as learning partners.
Gulbenkian Institute for Molecular Medicine (GIMM).

## Layout

| Folder | Contents |
|---|---|
| `notebooks/day1` … `day3` | Hands-on Colab notebooks (`DayN_HOk_topic.ipynb`) |
| `slides/day1` … `day3` | Interactive HTML slide decks (`DayN_Tk_topic.html`). Open in any browser, work offline |
| `slides/_build` | Sources and build script for the decks |
| `assets/gimm_logo` | Official GIMM logo files (use subject to GIMM approval) |
| `archive` | Superseded material kept for reference |

## Programme

- **Day 1 – Statistical inference foundations:** sampling variability, confidence intervals, p-values, false positives, power, choosing a test (GBSG2 breast cancer trial)
- **Day 2 – Regression:** tests as linear models, multiple regression, diagnostics, logistic and count models, differential expression and FDR (airway RNA-seq)
- **Day 3 – Survival analysis, experimental design, capstone** *(in preparation)*

## Presenting the slides

Open the `.html` file in a browser. ← → to navigate, **N** speaker notes, **F** fullscreen.
To rebuild after editing a source: `python slides/_build/build.py`.

## Facilitator notes

Kept in a separate private repository (expected results would spoil the in-class reveals). Facilitators clone it into `facilitator/` inside this folder; it is git-ignored here.
