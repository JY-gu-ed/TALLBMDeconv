# TALLBMDeconv

Reference-based paediatric T-ALL bone-marrow bulk RNA-seq deconvolution with optional audited flow-cytometry constraints and paired single-cell evaluation.

## v0.1.1 code-only prerelease

This repository distributes GPL-3 code, documentation and entirely synthetic examples. Patient data, clinical/FCS files and patient-derived single-cell references are excluded. Use a custom reference, or a separately authorised reference RDS. Bundled-reference loading functions report an actionable error when reference data are unavailable.

Features: reference construction/loading; RNA-only inference; optional flow observations; contribution and within-group expression extraction; QC; paired benchmark evaluation. RNA contribution, cell-scale proportions and within-group expression are distinct outputs. Gates supplied to a fit are calibration targets, not independent validation.

## Reproducible supported environment

Linux/amd64; R 4.5.2; Bioconductor 3.22; BayesPrism 2.2.3. Dockerfile pins the official Bioconductor image by digest. Public dependency sources, versions and SHA256 hashes are in `release/DEPENDENCY_SOURCE_LOCK.tsv` and `release/BAYESPRISM_SOURCE_LOCK.json`.

```sh
docker build --platform linux/amd64 --progress=plain -t tallbmdeconv:0.1.1 .
```

The build restores locked dependencies, installs the tested package archive, runs the synthetic custom-reference workflow and R CMD check. See GitHub Actions for the actual repository CI outcome. Container logs retain native diagnostic messages.

Alternatively, run `bash release/run_linux_apptainer.sh /new/empty/workdir`. The supported startup is `release/run_supported_example.R`, which explicitly loads BayesPrism at top level before the unchanged synthetic example. The example covers custom references, RNA-only, flow-enabled and disabled-flow fits, paired evaluation, extraction and QC. Short example fits test functionality, not statistical convergence.

## Verification and known issue

Prior clean-image Apptainer testing completed R CMD check with 0 errors, 0 warnings and 1 C++11 NOTE; testthat recorded 70 passes and one reference-data-dependent skip. Nested/lazy upstream namespace initialisation can emit native stack-imbalance diagnostics, including during testthat. These are retained, not suppressed. Top-level explicit loading avoided the observed diagnostics in the supported workflow; eight tested loading/GC variants produced identical numerical outputs. This is a startup workaround, not an upstream memory-safety fix.

The frozen 125-patient scientific deliverable and historical benchmark are separate controlled resources and are not rerun or distributed by this code release. Ref73 was the historical held-out benchmark/125 reference; Ref81 is a future deployment reference requiring separate distribution authorisation and does not inherit the independent-validation claim.

## Citation and licence

See `CITATION.cff`, `inst/CITATION`, `LICENSE` and `inst/THIRD_PARTY.md`. Cite BayesPrism and any reference data used. No DOI is claimed until an archive actually assigns one.
