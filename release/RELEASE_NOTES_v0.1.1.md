# TALLBMDeconv v0.1.1 — code-only prerelease

First public code release: custom single-cell reference construction/loading, bulk inference, optional audited flow observations, contribution/expression extraction, QC and paired evaluation.

Includes source code, the previously tested R package archive, synthetic examples, digest-pinned Linux environment and source/hash dependency locks. No patient data, clinical/FCS tables or patient-derived reference matrices are included. Users must supply their own or separately authorised reference.

Supported target: Linux/amd64, R 4.5.2, Bioconductor 3.22, BayesPrism 2.2.3. Release creation is gated on repository CI: locked environment build, synthetic workflow, R CMD check (zero errors/warnings required), and inference-source integrity checks. Full CI logs are retained as workflow artifacts.

Known issue: lazy/nested upstream namespace initialisation can emit native stack-imbalance diagnostics. The supported workflow explicitly loads BayesPrism at top level; diagnostics are not suppressed. Prior tests recorded identical numerical outputs across eight loading/GC variants, but this does not establish global memory safety. One authorised-reference test is intentionally skipped when those data are absent; a C++11 NOTE may remain.

This prerelease does not redistribute the atlas or the 125-patient research deliverable. A reference resource release and DOI archive are separate actions.
