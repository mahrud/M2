#!/bin/bash
#SBATCH --job-name=msolve-homog
#SBATCH --array=0-3
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=out/%x-%A_%a.out
#SBATCH --error=out/%x-%A_%a.err
#
# Homogeneous modules.  Sizes come from a calibration sweep and put msolve in
# the 3-15 second band, so no measurement is dominated by the fixed cost of
# handing the matrix to the engine and taking it back.
#
#   usage: mkdir -p out && sbatch --account=<your-account> bench-homog.sbatch
#   env:   M2  path to the Macaulay2 executable (default: M2 on PATH)

module load StdEnv/2023 gcc/12.3 boost/1.85.0 eigen/3.4.0 tbb/2021.10.0 2>/dev/null
cd "$SLURM_SUBMIT_DIR" || exit 1
M2=${M2:-M2}
echo "host=$(hostname) task=$SLURM_ARRAY_TASK_ID started=$(date -Is)"
case "$SLURM_ARRAY_TASK_ID" in
  0) ./drive.sh specs/hsat.txt    out/hsat.csv    400 time   ;;
  1) ./drive.sh specs/hmon.txt    out/hmon.csv    400 time   ;;
  2) ./drive.sh specs/hmel.txt    out/hmel.csv    400 time   ;;
  3) ./drive.sh specs/hverify.txt out/hverify.csv 500 verify ;;
esac
echo "task=$SLURM_ARRAY_TASK_ID finished=$(date -Is)"
