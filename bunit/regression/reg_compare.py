#!/usr/bin/env python3
# Metris: high-order metric-based non-manifold tetrahedral remesher
# Copyright (C) 2023-2025, Massachusetts Institute of Technology
# Licensed under The GNU Lesser General Public License, version 2.1
# See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php
"""Compare two sets of regression results (base, head).

Each set is a directory of <case>.json files, as written by the regression
cases under METRIS_REG_OUT (see reg_common.hxx). Prints a table with one row
per (case, run, metric) and a verdict: better, same, worse, new, missing,
failed, fixed, params.

Correctness (quality stats, run status, parameters) is one-sided: improvements
pass. Same tolerances as reg_common.hxx. Performance rows (--perf: cpu_s,
instructions) are reported, never part of the exit code.

Exit: 0 no correctness got worse, 1 correctness got worse, 2 bad input.
"""

import argparse
import json
import sys
from pathlib import Path

SCHEMA = 1  # keep in sync with REG_SCHEMA in reg_common.hxx
TOL_PCTUNIT = 1.0  # REG_TOL_PCTUNIT
TOL_QUALITY = 1.0e-2  # REG_TOL_QUALITY
# Relative change under which a perf metric is reported "same".
PERF_SAME = {"cpu_s": 0.05, "instructions": 0.005}

WORSE = {"worse", "missing", "failed", "params"}


def stat_check(key):
    """(tolerance, sign) if key is a checked stat, else None.

    sign +1: higher is better (% unit edges). -1: lower is better (quality:
    the stats hold the conformity error |q - 1|, 0 is ideal).
    """
    if key.startswith("pctunit"):
        return TOL_PCTUNIT, 1
    if key.startswith("avgquaP1_") or key.startswith("avgquaPk_"):
        return TOL_QUALITY, -1
    return None


class InputError(Exception):
    pass


def load_set(path):
    path = Path(path)
    if not path.is_dir():
        raise InputError(f"{path}: not a directory")
    cases = {}
    for f in sorted(path.glob("*.json")):
        try:
            data = json.loads(f.read_text())
        except (OSError, json.JSONDecodeError) as e:
            raise InputError(f"{f}: {e}")
        if not isinstance(data, dict) or data.get("schema") != SCHEMA:
            raise InputError(f"{f}: schema {data.get('schema') if isinstance(data, dict) else None}, expected {SCHEMA}")
        cases[data["case"]] = data
    if not cases:
        raise InputError(f"{path}: no result files")
    return cases


def one_value(cases, key, label):
    vals = {c.get(key) for c in cases.values()}
    if len(vals) != 1:
        raise InputError(f"{label}: several values of {key}: {sorted(map(str, vals))}")
    return vals.pop()


def compare(base, head, perf):
    rows = []

    def row(case, run, metric, b, h, verdict, kind="correctness"):
        rows.append(dict(case=case, run=run, metric=metric, base=b, head=h,
                         verdict=verdict, kind=kind))

    for case in sorted(set(base) | set(head)):
        if case not in head:
            row(case, "", "case", "present", None, "missing")
            continue
        if case not in base:
            row(case, "", "case", None, "present", "new")
            continue
        bruns, hruns = base[case]["runs"], head[case]["runs"]
        for run in sorted(set(bruns) | set(hruns)):
            if run not in hruns:
                row(case, run, "run", "present", None, "missing")
                continue
            if run not in bruns:
                row(case, run, "run", None, hruns[run]["status"], "new")
                continue
            b, h = bruns[run], hruns[run]
            if b["status"] != "ok" or h["status"] != "ok":
                if b["status"] == "ok":
                    verdict = "failed"
                elif h["status"] == "ok":
                    verdict = "fixed"
                else:
                    verdict = "same"
                row(case, run, "status", b["status"], h["status"], verdict)
                continue
            if b["params_nondefault"] != h["params_nondefault"]:
                diff = sorted(k for k in set(b["params_nondefault"]) | set(h["params_nondefault"])
                              if b["params_nondefault"].get(k) != h["params_nondefault"].get(k))
                row(case, run, "params", ",".join(diff), ",".join(diff), "params")
                continue
            bstats, hstats = b["stats"] or {}, h["stats"] or {}
            for key in sorted(bstats):
                check = stat_check(key)
                if check is None:
                    continue
                tol, sign = check
                if key not in hstats:
                    row(case, run, key, bstats[key], None, "missing")
                    continue
                bv, hv = bstats[key], hstats[key]
                gain = sign * (hv - bv)
                verdict = "worse" if gain < -tol else "better" if gain > tol else "same"
                row(case, run, key, bv, hv, verdict)
            if perf:
                for key, rel in PERF_SAME.items():
                    bv, hv = b.get(key), h.get(key)
                    if bv is None or hv is None or bv <= 0:
                        verdict = "n/a"
                    else:
                        ratio = hv / bv
                        verdict = "worse" if ratio > 1 + rel else "better" if ratio < 1 - rel else "same"
                    row(case, run, key, bv, hv, verdict, kind="perf")
    return rows


def fmt_val(v):
    if v is None:
        return "-"
    if isinstance(v, float):
        return f"{v:.4g}"
    return str(v)


def print_table(rows, out=sys.stdout):
    cols = ["case", "run", "metric", "base", "head", "ratio", "verdict"]
    lines = []
    for r in rows:
        ratio = ""
        if r["kind"] == "perf" and isinstance(r["base"], (int, float)) and isinstance(r["head"], (int, float)) and r["base"]:
            ratio = f"{r['head'] / r['base']:.3f}"
        lines.append([r["case"], r["run"], r["metric"], fmt_val(r["base"]), fmt_val(r["head"]), ratio, r["verdict"]])
    widths = [max([len(c)] + [len(l[i]) for l in lines]) for i, c in enumerate(cols)]
    print("  ".join(c.ljust(w) for c, w in zip(cols, widths)), file=out)
    for l in lines:
        print("  ".join(v.ljust(w) for v, w in zip(l, widths)), file=out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("base", help="directory of base result files")
    ap.add_argument("head", help="directory of head result files")
    ap.add_argument("--json", metavar="FILE", help="also write the rows and verdicts as JSON")
    ap.add_argument("--perf", action="store_true", help="add cpu_s and instructions rows")
    args = ap.parse_args(argv)

    try:
        base, head = load_set(args.base), load_set(args.head)
        meta = {}
        for key in ("hwid", "build_type"):
            bv, hv = one_value(base, key, args.base), one_value(head, key, args.head)
            if bv != hv:
                raise InputError(f"{key} differs: base {bv}, head {hv}")
            meta[key] = bv
        for name, s in (("base", base), ("head", head)):
            meta[name] = {"commit": one_value(s, "commit", name), "dirty": any(c.get("dirty") for c in s.values())}
    except InputError as e:
        print(f"reg_compare: {e}", file=sys.stderr)
        return 2

    rows = compare(base, head, args.perf)
    shown = rows if args.perf else [r for r in rows if r["kind"] == "correctness"]
    print_table(shown)
    worse = [r for r in rows if r["kind"] == "correctness" and r["verdict"] in WORSE]
    code = 1 if worse else 0
    print(f"\n{len(worse)} correctness row(s) worse." if worse else "\nNo correctness row worse.")

    if args.json:
        Path(args.json).write_text(json.dumps(dict(meta, exit=code, rows=rows), indent=2) + "\n")
    return code


if __name__ == "__main__":
    sys.exit(main())
