#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
work=${1:?Provide a NEW empty work directory}
mkdir -p "$work"
if [ -e "$work/R_library" ]; then echo "Use a new work directory" >&2; exit 1; fi
mkdir -p "$work/release" "$work/tmp"
cp "$root"/release/*.json "$root"/release/*.tsv "$root"/release/*.R "$work/release/"
cp "$root/dist/TALLBMDeconv_0.1.1.tar.gz" "$root/examples/PUBLIC_EXAMPLE.R" "$work/"
image='docker://bioconductor/bioconductor_docker@sha256:67ab741591c5cd068f8993ca5f0371794c300da9cc5e0b3b2c1c90c9058e4c6f'
apptainer exec --cleanenv --containall --no-home --no-mount cwd,hostfs --bind "$work:/work,$work/tmp:/tmp" --env TMPDIR=/work/tmp --pwd /work "$image" Rscript --vanilla /work/release/restore_locked.R
apptainer exec --cleanenv --containall --no-home --no-mount cwd,hostfs --bind "$work:/work,$work/tmp:/tmp" --env TMPDIR=/work/tmp --pwd /work "$image" Rscript --vanilla /work/release/run_supported_example.R synthetic_run
