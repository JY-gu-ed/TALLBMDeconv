# TALLBMDeconv 0.1.1

* Verify modeled counts and effective flow inputs before reusing benchmark fits;
  missing supplied fits are not silently recomputed.
* Handle partial/empty per-patient flow tables. Strength-zero observations are
  recorded as requested but unused.
* Fix and record the frozen RNG kind and restore caller RNG state.
* Reject parent/component namespace collisions, invalid truth mappings and
  incomplete pseudobulk truth. Export truth coverage and propagate file errors.
* Mark completion only after QA; reject incomplete result reads.
* Frozen references, numerical sampler kernels and 125-patient matrices unchanged.

# TALLBMDeconv 0.1.0

* First functional package: default/custom reference, two-chain fixed-reference
  initial allocation, typed optional flow, file-backed source/expression output,
  same-fit conservation QA and patient-level paired benchmark.
* Distinct benchmark73 and fullAtlas81 resources in the authorised distribution.
* Standardised TALL125_DECONV_V1 released separately without refitting125.
* Public reference redistribution remains subject to author/source approval.
