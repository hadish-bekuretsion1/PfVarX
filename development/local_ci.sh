#!/usr/bin/env bash
# Local CI: the same checks as .github/workflows/ci.yaml, run on this computer before every
# merge and release (GitHub CI is started by hand only). Development material: never released.
#
# Usage, from anywhere:
#   development/local_ci.sh            full run, about 17 minutes
#   development/local_ci.sh --quick    skip the integration tests, about 3 minutes
# The log goes to $PFVARX_CI_LOG_DIR (default ~/pfvarx-lab) as local_ci_<commit>.log.
# Run one copy at a time: two runs of the same commit write the same log.
set -euo pipefail
cd "$(dirname "$0")/.."
logdir=${PFVARX_CI_LOG_DIR:-$HOME/pfvarx-lab}
mkdir -p "$logdir"
log=$logdir/local_ci_$(git rev-parse --short HEAD).log
step() { printf '\n== %s\n' "$1"; }
{
  echo "Branch $(git branch --show-current) at $(git rev-parse --short HEAD), $(date)"
  if [ -n "$(git status --porcelain)" ]; then echo "Note: uncommitted changes are included in this run"; fi
  step "Environment from the lock file";  pixi install --locked
  step "Lint";                            pixi run lint
  step "Type checking";                   pixi run typecheck
  step "Import and help smoke";           pixi run pfvarx --help > /dev/null && pixi run pfvarx version
  step "Unit and scientific tests";       pixi run pytest -q tests/unit tests/scientific
  if [ "${1:-}" != "--quick" ]; then
    step "Integration tests";             pixi run pytest -q tests/integration
  fi
  step "Software-only source archive"
  out=$(mktemp -d)
  version="$(pixi run --locked pfvarx version)"
  pixi run python scripts/package.py --output "$out" > /dev/null
  (cd "$out" && sha256sum -c "pfvarx-${version}-source.zip.sha256" && python3 -m zipfile -e "pfvarx-${version}-source.zip" unpacked && python3 unpacked/pfvarx/scripts/package.py --check)
  rm -rf "$out"
  echo; echo "LOCAL CI PASSED"
} 2>&1 | tee "$log"
