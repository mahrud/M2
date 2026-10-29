#!/bin/bash
# Run every measurement in a spec file, one Macaulay2 process each, appending a
# CSV line per measurement.
#
# A separate process per measurement is what keeps a Groebner basis cached by
# one method out of the next method's timing, and it lets the timeout be
# enforced from outside: msolve does not poll for interrupts inside its own C
# loops, so an in-process alarm cannot stop it.
#
#   usage: drive.sh <specfile> <outcsv> <per-case-timeout-seconds> [time|verify]
#   env:   M2   path to the Macaulay2 executable (default: M2 on PATH)
set -u
SPEC=$1; OUT=$2; TMO=$3; MODE=${4:-time}
M2=${M2:-M2}
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
SCRIPT=runone.m2; PREFIX=RESULT
if [ "$MODE" = verify ]; then SCRIPT=verify.m2; PREFIX=VERIFY; fi
mkdir -p "$(dirname "$OUT")"
: > "$OUT"
while read -r fam cs meth thr; do
  case "${fam:-}" in ''|\#*) continue;; esac
  line=$(FAMILY="$fam" CASE="$cs" METHOD="${meth:-}" THREADS="${thr:-1}" MODE="$MODE" \
         timeout -s KILL "$TMO" "$M2" --script "$SCRIPT" 2>/dev/null | grep "^$PREFIX," | head -1)
  if [ -z "$line" ]; then
    line="$PREFIX,$fam,$cs,${meth:-},${thr:-1},$TMO,$TMO,timeout,killed-at-${TMO}s"
  fi
  printf '%s\n' "$line" | tee -a "$OUT"
done < "$SPEC"
