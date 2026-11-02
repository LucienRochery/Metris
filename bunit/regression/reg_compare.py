#!/usr/bin/env python3
# Metris: high-order metric-based non-manifold tetrahedral remesher
# Copyright (C) 2023-2025, Massachusetts Institute of Technology
# Licensed under The GNU Lesser General Public License, version 2.1
# See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php
"""Compare two sets of regression results (base, head).

Each set is a directory of <case>.json files, as written by the regression
cases under METRIS_REG_OUT (see reg_common.hxx). Prints a table with one row
per (case, run, metric) and a verdict: better, same, worse, new, missing,
failed, fixed, params. A new run that throws is "failed".

Correctness (quality stats, run status, parameters) is one-sided: improvements
pass. Same tolerances as reg_common.hxx. Performance rows (--perf: cpu_s,
instructions) are reported, never part of the exit code.

Exit: 0 no correctness got worse, 1 correctness got worse, 2 bad input
(--json is then written with exit 2 and an error message).
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


def _reject_constant(name):
    raise ValueError(f"non-finite number {name}")


def _check_number(v, where):
    if v is not None and (isinstance(v, bool) or not isinstance(v, (int, float))):
        raise InputError(f"{where}: expected a number or null, got {type(v).__name__}")


def validate(data, f):
    """Raise InputError unless data looks like a result file."""
    if not isinstance(data, dict) or data.get("schema") != SCHEMA:
        got = data.get("schema") if isinstance(data, dict) else None
        raise InputError(f"{f}: schema {got}, expected {SCHEMA}")
    if not isinstance(data.get("case"), str):
        raise InputError(f"{f}: no case name")
    if data["case"] != f.stem:
        raise InputError(f"{f}: holds case {data['case']}, not {f.stem}")
    runs = data.get("runs")
    if not isinstance(runs, dict):
        raise InputError(f"{f}: runs is not an object")
    for name, run in runs.items():
        where = f"{f}: run {name}"
        if not isinstance(run, dict) or not isinstance(run.get("status"), str):
            raise InputError(f"{where}: no status")
        if run["status"] != "ok":
            continue
        if not isinstance(run.get("params_nondefault"), dict):
            raise InputError(f"{where}: params_nondefault is not an object")
        if not isinstance(run.get("stats"), dict):
            raise InputError(f"{where}: stats is not an object")
        for key, v in run["stats"].items():
            _check_number(v, f"{where}: stat {key}")
        for key in PERF_SAME:
            _check_number(run.get(key), f"{where}: {key}")


def load_set(path):
    path = Path(path)
    if not path.is_dir():
        raise InputError(f"{path}: not a directory")
    cases = {}
    for f in sorted(path.glob("*.json")):
        try:
            data = json.loads(f.read_text(encoding="utf-8"), parse_constant=_reject_constant)
        except (OSError, UnicodeDecodeError, ValueError) as e:
            raise InputError(f"{f}: {e}")
        validate(data, f)
        # File name = case name (checked above), so no two files share a case.
        cases[data["case"]] = data
    if not cases:
        raise InputError(f"{path}: no result files")
    return cases


def one_value(cases, key, label):
    vals = {c.get(key) for c in cases.values()}
    if len(vals) != 1:
        raise InputError(f"{label}: several values of {key}: {sorted(map(str, vals))}")
    return vals.pop()


def set_info(cases):
    """Commit(s) and dirty flag of a result set. Files may come from different
    commits (e.g. baselines, updated case by case): commit is then null."""
    commits = sorted({str(c.get("commit")) for c in cases.values()})
    return {"commit": commits[0] if len(commits) == 1 else None,
            "commits": commits,
            "dirty": any(c.get("dirty") for c in cases.values())}


def compare(base, head, perf):
    rows = []

    def row(case, run, metric, b, h, verdict, kind="correctness"):
        rows.append(dict(case=case, run=run, metric=metric, base=b, head=h,
                         verdict=verdict, kind=kind))

    for case in sorted(set(base) | set(head)):
        if case not in head:
            row(case, "", "case", "present", None, "missing")
            continue
        bruns = base[case]["runs"] if case in base else {}
        hruns = head[case]["runs"]
        if case not in base:
            row(case, "", "case", None, "present", "new")
        for run in sorted(set(bruns) | set(hruns)):
            if run not in hruns:
                row(case, run, "run", "present", None, "missing")
                continue
            h = hruns[run]
            if run not in bruns:
                # A new run that throws fails, as in reg_common.hxx.
                row(case, run, "run", None, h["status"], "new" if h["status"] == "ok" else "failed")
                continue
            b = bruns[run]
            if b["status"] != "ok" or h["status"] != "ok":
                if b["status"] == "ok":
                    verdict = "failed"
                elif h["status"] == "ok":
                    verdict = "fixed"
                else:
                    verdict = "same"
                row(case, run, "status", b["status"], h["status"], verdict)
                continue
            bp, hp = b["params_nondefault"], h["params_nondefault"]
            # Typed comparison, as nlohmann::json's (true != 1 there, not in Python).
            def same(x, y):
                return json.dumps(x, sort_keys=True) == json.dumps(y, sort_keys=True)
            if not same(bp, hp):
                diff = sorted(k for k in set(bp) | set(hp) if k not in bp or k not in hp or not same(bp[k], hp[k]))
                row(case, run, "params", {k: bp.get(k) for k in diff},
                    {k: hp.get(k) for k in diff}, "params")
                continue
            bstats, hstats = b["stats"], h["stats"]
            for key in sorted(bstats):
                check = stat_check(key)
                if check is None:
                    continue
                tol, sign = check
                bv, hv = bstats[key], hstats.get(key)
                # null = NaN on the C++ side (e.g. % unit edges with no edges).
                if hv is None:
                    verdict = "same" if bv is None and key in hstats else "missing"
                elif bv is None:
                    verdict = "new"
                else:
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
    if isinstance(v, dict):
        return ",".join(f"{k}={json.dumps(x)}" for k, x in v.items())
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


def run(args):
    """Returns (exit code, JSON report)."""
    try:
        base, head = load_set(args.base), load_set(args.head)
        meta = {}
        for key in ("hwid", "build_type"):
            bv, hv = one_value(base, key, args.base), one_value(head, key, args.head)
            if bv != hv:
                raise InputError(f"{key} differs: base {bv}, head {hv}")
            meta[key] = bv
    except InputError as e:
        print(f"reg_compare: {e}", file=sys.stderr)
        return 2, {"exit": 2, "error": str(e), "rows": []}
    meta["base"], meta["head"] = set_info(base), set_info(head)

    rows = compare(base, head, args.perf)
    shown = rows if args.perf else [r for r in rows if r["kind"] == "correctness"]
    print_table(shown)
    worse = [r for r in rows if r["kind"] == "correctness" and r["verdict"] in WORSE]
    code = 1 if worse else 0
    print(f"\n{len(worse)} correctness row(s) worse." if worse else "\nNo correctness row worse.")
    return code, dict(meta, exit=code, rows=rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("base", help="directory of base result files")
    ap.add_argument("head", help="directory of head result files")
    ap.add_argument("--json", metavar="FILE", help="also write the rows and verdicts as JSON (on every exit)")
    ap.add_argument("--perf", action="store_true", help="add cpu_s and instructions rows")
    args = ap.parse_args(argv)

    try:
        code, report = run(args)
    except Exception as e:  # anything unforeseen is bad input, never "worse"
        print(f"reg_compare: internal error: {type(e).__name__}: {e}", file=sys.stderr)
        code, report = 2, {"exit": 2, "error": f"{type(e).__name__}: {e}", "rows": []}

    if args.json:
        try:
            Path(args.json).write_text(json.dumps(report, indent=2) + "\n")
        except OSError as e:
            print(f"reg_compare: cannot write {args.json}: {e}", file=sys.stderr)
            return 2
    return code


if __name__ == "__main__":
    sys.exit(main())
