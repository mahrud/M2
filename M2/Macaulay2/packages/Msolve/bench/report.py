#!/usr/bin/env python3
"""Turn the per-measurement CSV lines into per-case comparisons."""
import csv, sys, os
from collections import defaultdict

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

def load(fn):
    rows = []
    p = os.path.join(OUT, fn)
    if not os.path.exists(p): return rows
    for line in open(p):
        parts = line.strip().split(",")
        if len(parts) < 9 or parts[0] != "RESULT": continue
        _, fam, case, meth, thr, secs, cpu, status, detail = parts[:9]
        try: secs = float(secs)
        except ValueError: secs = float("nan")
        rows.append(dict(fam=fam, case=case, meth=meth, secs=secs,
                         status=status, detail=detail))
    return rows

def fmt(x):
    if x != x: return "   n/a"
    if x >= 100: return f"{x:6.0f}"
    if x >= 10:  return f"{x:6.1f}"
    if x >= 1:   return f"{x:6.2f}"
    return f"{x:6.3f}"

def table(rows, methods, title, baseline, contender):
    if not rows: return
    by = defaultdict(dict)
    stat = defaultdict(dict)
    order = []
    for r in rows:
        if r["case"] not in by: order.append(r["case"])
        by[r["case"]][r["meth"]] = r["secs"]
        stat[r["case"]][r["meth"]] = r["status"]
    w = max(len(c) for c in order) + 1
    print(f"\n{title}")
    print("-" * (w + 9 * len(methods) + 12))
    print("case".ljust(w) + "".join(m.rjust(9) for m in methods) + "  speedup")
    print("-" * (w + 9 * len(methods) + 12))
    speedups = []
    for c in order:
        line = c.ljust(w)
        for m in methods:
            v = by[c].get(m, float("nan"))
            s = stat[c].get(m, "")
            line += (">" + fmt(v).strip()).rjust(9) if s == "timeout" else fmt(v).rjust(9)
        b, k = by[c].get(baseline), by[c].get(contender)
        if b and k and k == k and b == b and k > 0:
            sp = b / k
            tb = stat[c].get(baseline) == "timeout"
            speedups.append(sp)
            line += f"  {sp:7.1f}x" + ("+" if tb else "")
        else:
            line += "      n/a"
        print(line)
    if speedups:
        speedups.sort()
        med = speedups[len(speedups)//2]
        print("-" * (w + 9 * len(methods) + 12))
        print(f"median speedup ({baseline} / {contender}): {med:.1f}x"
              f"   min {min(speedups):.1f}x   max {max(speedups):.1f}x")

sat = load("sat.csv")
table(sat, ["msolve", "msolve-guarded", "m2-iterate", "m2-default"],
      "MODULE SATURATION   seconds, single thread ('>' = hit the timeout)",
      "m2-iterate", "msolve")
table(sat, ["msolve", "msolve-guarded", "m2-iterate", "m2-default"],
      "MODULE SATURATION, guarded msolve vs Macaulay2", "m2-iterate", "msolve-guarded")

elim = load("elim.csv")
table(elim, ["msolve", "msolve-exp", "m2"],
      "IDEAL ELIMINATION   seconds, single thread", "m2", "msolve")

mel = load("melim.csv")
table(mel, ["msolve", "m2"],
      "MODULE ELIMINATION   seconds, single thread", "m2", "msolve")

vp = os.path.join(OUT, "verify.csv")
if os.path.exists(vp):
    agree, differ, timeout, declined = 0, [], [], []
    for line in open(vp):
        p = line.strip().split(",")
        if len(p) < 4 or p[0] != "VERIFY": continue
        if "timeout" in line:      timeout.append(f"{p[1]}/{p[2]}")
        elif p[3] == "agree":
            agree += 1
            if len(p) > 4 and "declined" in p[4]: declined.append(f"{p[1]}/{p[2]}")
        else: differ.append(line.strip())
    print(f"\nVERIFICATION")
    print("-" * 56)
    print(f"  agree with Macaulay2 : {agree}")
    print(f"  disagree             : {len(differ)}")
    print(f"  not decided (timeout): {len(timeout)}" +
          ("  [" + ", ".join(timeout) + "]" if timeout else ""))
    if declined: print("  msolve declined, fell back: " + ", ".join(declined))
    for b in differ: print("  DISAGREEMENT: " + b)
