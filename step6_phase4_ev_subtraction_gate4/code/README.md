# Pillar 4 — code 

Python code for Pillar 4

## Layout

* `synthesis/` — the Pillar-4 subtraction + gate proper.
* `synthesis\_tables\_figures/` — the cross-pillar table/figure builders fed by a
single data layer. (Two folders avoid a name collision; both write to
`04\_synthesis/data/` and `04\_synthesis/figures/`.)

## Run order

1. `synthesis/pillar4\_synthesis.py` — whole-landscape EV subtraction: route
composition of the 34 reproducible survivors (24 ECM + 10 soluble +
0 EV-compatible) and the 3 cross-atlas consensus axes.
2. `synthesis/gate4\_checks.py` — Gate 4 synthesis-integrity checks (internal
soundness + faithfulness to frozen artifacts and locked Gate 1–3 verdicts).
3. `synthesis/pillar4\_figures.py` — Pillar-4 figures (F\_FUNNEL evidence
refinement, F6 subtraction, F7 route map, cross-pillar panels).
4. `synthesis\_tables\_figures/synthesis\_data.py` — data layer: `load\_all()`
parses every frozen artifact; `verify()` asserts the re-derivation
identities. Run standalone to audit.
5. `synthesis\_tables\_figures/build\_tables.py` — synthesis tables T1–T5.
6. `synthesis\_tables\_figures/figures\_tier1.py` + `figures\_tier2.py` —
headline cross-pillar synthesis figures (A–D and E–H).

## Notes

* Absolute paths recorded inside scripts reflect the original analysis
environment; adapt them before re-running.

