# Figure QA summary (Development figures, 180 mm wide)

| Figure | Alignment gate | Collision audit | Glyph floor |
|---|---|---|---|
| Fig1 | PASS | PASS | see qa/Fig1.pdf-text-audit.json |
| Fig2 | PASS | PASS | see qa/Fig2.pdf-text-audit.json |
| Fig3 | PASS | PASS | see qa/Fig3.pdf-text-audit.json |
| Fig4 | PASS | PASS | see qa/Fig4.pdf-text-audit.json |
| FigS1 | PASS | PASS | see qa/FigS1.pdf-text-audit.json |

Source preflight (validate_figure.py): READY FOR VISUAL QA for all scripts. Remaining WARNs: DATA-EXCLUSION (missing-value drops are counted and printed by each script), DEMO-DATA (random jitter used only for strip-plot display), LOG-GUARD (pseudocount of 1 applied before log2).
Visual audit at final size: each panel inspected; no overlaps; all text 8 pt or larger except the log2 subscript (5.6 pt).
