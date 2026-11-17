#!/usr/bin/env bash
set -euo pipefail
cd /work
mkdir -p ci_qa
export R_LIBS_USER=/work/R_library
export TMPDIR=/work/tmp
R CMD check --no-manual TALLBMDeconv_0.1.1.tar.gz >ci_qa/R_CMD_check.stdout 2>ci_qa/R_CMD_check.stderr
cp -r TALLBMDeconv.Rcheck ci_qa/
if grep -Eq '[1-9][0-9]* (ERROR|WARNING)' TALLBMDeconv.Rcheck/00check.log; then
  cat TALLBMDeconv.Rcheck/00check.log
  exit 1
fi
cp RESTORE_ACCEPTANCE.txt ci_qa/
printf 'PASS_REPOSITORY_CONTAINER_CHECK
' >ci_qa/ACCEPTANCE.txt
