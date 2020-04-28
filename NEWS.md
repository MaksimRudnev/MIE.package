
# Large update to 0.5-1 ------------------

*April, 28, 2020*

Several new functions are introduced or moved from other packages/gists.

- automate running measurement invariance alignment in Mplus; extract its results; run simulations on alignment and extract its key results back to R objects;
- `measurementInvarianceMplus` creates data, input, runs it in Mplus, and extracts its results back to R; now handles categorical variables as well. Still at early stage.
- `globalMI` now can handle categorical variables, applies proper methods to them; still at early stage.
- `mgcfa_diagnose` prints a quick summary of a multiple group CFA models fitted with lavaan; it includes key fit indices, aggregated across modification indices, and (optionally) score test of constraints (lavTestScore with groups labels),