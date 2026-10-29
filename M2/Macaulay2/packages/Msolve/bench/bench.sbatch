#!/bin/bash
#SBATCH --job-name=msolve-bench
#SBATCH --array=0-3
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=out/%x-%A_%a.out
#SBATCH --error=out/%x-%A_%a.err
#
# One array task per benchmark family.  Single core throughout: Macaulay2's
# Groebner engine is serial, so msolve is held to one thread as well and the
# comparison is like for like.
#
#   usage: mkdir -p out && sbatch --account=<your-account> bench.sbatch
#   env:   M2  path to the Macaulay2 executable (default: M2 on PATH)

module load StdEnv/2023 gcc/12.3 boost/1.85.0 eigen/3.4.0 tbb/2021.10.0 2>/dev/null

cd "$SLURM_SUBMIT_DIR" || exit 1
M2=${M2:-M2}
MISSING=$(ldd "$(command -v "$M2" || echo "$M2")" 2>/dev/null | grep 'not found' || true)
if [ -n "$MISSING" ]; then echo "ABORT: missing libraries on compute node:"; echo "$MISSING"; exit 1; fi

echo "host=$(hostname) task=$SLURM_ARRAY_TASK_ID started=$(date -Is)"

case "$SLURM_ARRAY_TASK_ID" in
  0) ./drive.sh specs/sat.txt    out/sat.csv    120 time   ;;
  1) ./drive.sh specs/elim.txt   out/elim.csv   120 time   ;;
  2) ./drive.sh specs/melim.txt  out/melim.csv  120 time   ;;
  3) ./drive.sh specs/verify.txt out/verify.csv 240 verify ;;
esac

echo "task=$SLURM_ARRAY_TASK_ID finished=$(date -Is)"
