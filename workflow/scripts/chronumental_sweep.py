#!/usr/bin/env python3
"""Sweep chronumental's date-prior settings and measure what each one does (#146).

#145 confirmed that chronumental buys leaf-date accuracy with inflated branch
durations on both v11 production trees: the fitted rate lands 4.3-4.4x below
chronumental's own pre-fit root-to-tip estimate, ``length_cor`` falls monotonically
from 0.995 to ~0.38, and 57% of branches carry zero mutations while holding 46% of
the tree's total elapsed time. This script runs the experiment that decides what to
do about it.

Three design choices are load-bearing.

**It runs the shipped QC, not a copy of its metrics.** Every per-arm number except
the stdout parse comes from executing the real ``scripts/chronumental_qc.py`` against
that arm's ``node_dates.tsv`` -- including the tip-residual gate verdict, which is
the shipped ``record()``/``SystemExit`` path rather than a threshold re-applied here.
A sweep that re-implemented "median absolute tip residual" could disagree with the
gate the pipeline actually enforces, and would then be measuring the wrong quantity
with total confidence. The QC script writes its report *before* raising, so a
gate-failing arm still yields a full row -- breaching the gate is a result, not an
error.

**It never touches the source tree.** Both chronumental inputs
(``mat/filtered_tree.nwk``, ``mat/filtered_metadata.tsv.gz``) already exist in the
output tree, so each arm reads them and writes to its own scratch directory. No
pipeline rule runs, no production output is overwritten, and the baseline arm's job
is to reproduce numbers that are already on disk -- which is what validates the
harness before anything is concluded from the other arms.

**The baseline arm is run, not assumed.** ``A0`` passes no prior flags at all, so it
must reproduce the shipped fit. If it does not, the harness is wrong and every other
arm is uninterpretable.

Usage::

    # one arm at a time, or several concurrently (each ~22 GB / ~2 h on public v11)
    python workflow/scripts/chronumental_sweep.py run \\
        --tree-dir results/public_2026-08-06_v11 --out-dir <scratch> --jobs 6

    # then, serially (each ~40 GB): the real chronumental_qc per arm
    python workflow/scripts/chronumental_sweep.py qc \\
        --tree-dir results/public_2026-08-06_v11 --out-dir <scratch>

    # finally, the comparison table
    python workflow/scripts/chronumental_sweep.py report --out-dir <scratch>
"""

import argparse
import json
import re
import subprocess
import sys
import time
from importlib import util as importlib_util
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

#: The six arms of #146's deliverable 2. ``--steps`` is held at the pipeline's
#: configured 1000 for every arm so each is directly comparable to the shipped fit;
#: step sensitivity is deliberately a separate, single-arm follow-up (D4), since
#: ``length_cor`` is still falling at step 999 and confounding the two would make a
#: prior effect and an under-convergence effect indistinguishable.
ARMS = {
    "A0_baseline": [],
    "A1_variance_dates_3": ["--variance_dates", "3"],
    "A2_variance_dates_10": ["--variance_dates", "10"],
    "A3_variance_dates_30": ["--variance_dates", "30"],
    "B0_horseshoe_tau3.2": ["--model", "HorseShoeLike"],
    "B1_horseshoe_tau32": ["--model", "HorseShoeLike", "--initial_tau", "32"],
}

#: Windowing parameters the durations feed into. Defaults match config/config.yaml;
#: the fraction of branches longer than these is the direct downstream symptom.
WINDOW_SIZE_DAYS = 90
WINDOW_STEP_DAYS = 45


def load_module(name):
    """Import a workflow module by path, the way the pipeline scripts do."""
    spec = importlib_util.spec_from_file_location(
        name, REPO_ROOT / "workflow" / "modules" / f"{name}.py"
    )
    module = importlib_util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def conda_prefix_with(*packages):
    """The Snakemake-built conda env providing every one of ``packages``.

    The sweep runs outside Snakemake, so it has to find the same interpreters the
    pipeline would use. Resolved by looking rather than hardcoded: env hashes change
    whenever an env file does. Every requirement is checked, not just the headline
    one -- ``bte`` ships as an extension module (``bte.cpython-*.so``) and lives in
    more than one env here, only some of which also carry pandas, so matching on the
    first name alone picks an env the QC script then dies in.
    """
    for prefix in sorted((REPO_ROOT / ".snakemake" / "conda").glob("*_")):
        if all(
            list(prefix.glob(f"lib/python*/site-packages/{package}*"))
            for package in packages
        ):
            return prefix
    raise SystemExit(
        f"no Snakemake conda env provides all of {packages}; run the pipeline once "
        f"with --use-conda first so the envs exist"
    )


# --- run: one chronumental fit per arm ---------------------------------------


def run_arms(tree_dir, out_dir, arms, steps, jobs):
    tree = tree_dir / "mat" / "filtered_tree.nwk"
    dates = tree_dir / "mat" / "filtered_metadata.tsv.gz"
    for path in (tree, dates):
        if not path.exists():
            raise SystemExit(f"missing chronumental input: {path}")
    chronumental = conda_prefix_with("chronumental") / "bin" / "chronumental"

    running = []
    for arm in arms:
        arm_dir = out_dir / arm
        arm_dir.mkdir(parents=True, exist_ok=True)
        if (arm_dir / "node_dates.tsv").exists():
            print(f"[{arm}] already fitted, skipping")
            continue
        cmd = [
            str(chronumental),
            "--tree", str(tree),
            "--dates", str(dates),
            "--steps", str(steps),
            *ARMS[arm],
            "--dates_out", str(arm_dir / "node_dates.tsv"),
            "--tree_out", str(arm_dir / "timetree.nwk"),
        ]
        (arm_dir / "cmd.txt").write_text(" ".join(cmd) + "\n")
        stdout = open(arm_dir / "chronumental.stdout", "w")
        stderr = open(arm_dir / "chronumental.stderr", "w")
        # /usr/bin/time -v gives peak RSS and wall-clock without wrapping the fit in
        # a Python resource accounting layer that would only measure this process.
        proc = subprocess.Popen(
            ["/usr/bin/time", "-v", *cmd], stdout=stdout, stderr=stderr
        )
        print(f"[{arm}] pid {proc.pid}: {' '.join(ARMS[arm]) or '(no prior flags)'}")
        running.append((arm, proc, stdout, stderr, time.time()))
        while len([p for _, p, _, _, _ in running if p.poll() is None]) >= jobs:
            time.sleep(30)

    for arm, proc, stdout, stderr, started in running:
        code = proc.wait()
        stdout.close()
        stderr.close()
        elapsed = time.time() - started
        (out_dir / arm / "wall_seconds.txt").write_text(f"{elapsed:.1f}\n")
        print(f"[{arm}] exit {code} after {elapsed / 3600:.2f} h")


# --- qc: the real chronumental_qc per arm ------------------------------------


def qc_arms(tree_dir, out_dir, arms, residual_threshold):
    tree = tree_dir / "mat" / (
        "filtered_tree.pb"
        if (tree_dir / "mat" / "filtered_tree.pb").exists()
        else "filtered_tree_raw.pb"
    )
    metadata = tree_dir / "mat" / "filtered_metadata.tsv.gz"
    script = REPO_ROOT / "workflow" / "scripts" / "chronumental_qc.py"
    python = conda_prefix_with("bte", "pandas", "numpy") / "bin" / "python"

    for arm in arms:
        arm_dir = out_dir / arm
        if not (arm_dir / "node_dates.tsv").exists():
            print(f"[{arm}] no node_dates.tsv yet, skipping qc")
            continue
        qc_dir = arm_dir / "qc"
        qc_dir.mkdir(parents=True, exist_ok=True)
        if (qc_dir / "node_dates_qc.tsv").exists():
            print(f"[{arm}] qc already done, skipping")
            continue
        # A subprocess, not an in-process exec: the QC script needs the bte env's
        # interpreter, and each arm's ~40 GB of arrays should be reclaimed by the OS
        # between arms rather than accumulating in one process.
        driver = qc_dir / "_run_qc.py"
        driver.write_text(
            "import runpy, sys\n"
            f"sys.path.insert(0, {str(script.parent)!r})\n"
            "class S:\n"
            "    pass\n"
            "class NS(dict):\n"
            "    __getattr__ = dict.__getitem__\n"
            "s = S()\n"
            f"s.input = NS(tree={str(tree)!r}, node_dates={str(arm_dir / 'node_dates.tsv')!r}, "
            f"metadata={str(metadata)!r})\n"
            f"s.output = NS(report={str(qc_dir / 'node_dates_qc.tsv')!r}, "
            f"residuals={str(qc_dir / 'node_dates_residuals.csv')!r}, "
            f"branch_dt={str(qc_dir / 'node_dates_branch_dt.csv')!r}, "
            f"branch_dt_summary={str(qc_dir / 'node_dates_branch_dt_summary.csv')!r}, "
            f"internal_checks={str(qc_dir / 'node_dates_internal_checks.csv')!r})\n"
            f"s.params = NS(max_median_abs_tip_residual_days={residual_threshold!r}, "
            f"chronumental_stdout={str(arm_dir / 'chronumental.stdout')!r})\n"
            f"s.log = [{str(qc_dir / 'chronumental_qc.log')!r}]\n"
            "import builtins\n"
            "builtins.snakemake = s\n"
            f"runpy.run_path({str(script)!r}, init_globals={{'snakemake': s}}, "
            "run_name='__main__')\n"
        )
        print(f"[{arm}] running chronumental_qc...")
        started = time.time()
        proc = subprocess.run(
            [str(python), str(driver)],
            capture_output=True,
            text=True,
        )
        # A non-zero exit is expected when an arm breaches the tip-residual gate:
        # the report is written before the SystemExit, so the row is still complete
        # and the gate verdict is the shipped code's, not ours.
        (qc_dir / "qc_exit.json").write_text(
            json.dumps(
                {
                    "returncode": proc.returncode,
                    "stderr_tail": proc.stderr[-2000:],
                    "seconds": round(time.time() - started, 1),
                },
                indent=2,
            )
        )
        print(f"[{arm}] qc exit {proc.returncode} in {(time.time() - started) / 60:.1f} min")


# --- report: one row per arm -------------------------------------------------


def parse_stdout(path, node_date_qc):
    """Root-to-tip rate, the final step line, and whether length_cor was falling."""
    if not path.exists():
        return {}
    text = path.read_text()
    steps = re.findall(r"^step:(\d+)\t(.*)$", text, flags=re.MULTILINE)
    out = {"root_to_tip_rate": node_date_qc.parse_root_to_tip_rate(str(path))}
    if not steps:
        return out
    fields = dict(
        kv.split(":", 1) for kv in steps[-1][1].split("\t") if ":" in kv
    )
    out.update(
        {
            "final_step": int(steps[-1][0]),
            "fitted_rate": float(fields.get("mutation_rate", "nan")),
            "length_cor": float(fields.get("length_cor", "nan")),
            "date_error_med": float(fields.get("date_error_med", "nan")),
            "max_date_error": float(fields.get("max_date_error", "nan")),
        }
    )
    cors = [float(dict(kv.split(":", 1) for kv in s[1].split("\t") if ":" in kv)
                  ["length_cor"]) for s in steps]
    out["length_cor_first"] = cors[0]
    out["length_cor_falling_intervals"] = sum(
        1 for a, b in zip(cors, cors[1:]) if b < a
    )
    out["length_cor_intervals"] = len(cors) - 1
    return out


def parse_qc(qc_dir):
    """Pull the checks that decide the outcome out of the shipped report."""
    report = qc_dir / "node_dates_qc.tsv"
    if not report.exists():
        return {}
    rows = {}
    for line in report.read_text().splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) >= 3:
            rows[parts[0]] = (parts[1], parts[2])
    out = {}
    if "tip_residual" in rows:
        status, detail = rows["tip_residual"]
        out["tip_residual_gate"] = status
        for key, pattern in (
            ("tip_residual_median_d", r"median\|res\|=([\d.]+)d"),
            ("tip_residual_p95_d", r"p95\|res\|=([\d.]+)d"),
        ):
            m = re.search(pattern, detail)
            if m:
                out[key] = float(m.group(1))
    if "root_date_plausible" in rows:
        detail = rows["root_date_plausible"][1]
        m = re.search(r"dated (\d{4}-\d{2}-\d{2})", detail)
        out["root_date"] = m.group(1) if m else None
        out["root_date_inside_window"] = "inside" in detail
    if "branch_dt_vs_mutations" in rows:
        detail = rows["branch_dt_vs_mutations"][1]
        for key, pattern in (
            ("implied_subs_per_year", r"= ([\d.]+) subs/yr"),
            ("zero_mut_branch_time_pct", r"hold ([\d.]+)% of total branch-time"),
            ("incremental_slope_d_per_mut", r"slope (\d+) d/mutation"),
        ):
            m = re.search(pattern, detail)
            if m:
                out[key] = float(m.group(1))
    return out


def window_overruns(qc_dir):
    """Fraction of branches longer than one window, and than one step."""
    path = qc_dir / "node_dates_branch_dt.csv"
    if not path.exists():
        return {}
    import csv

    total = over_window = over_step = 0
    with open(path) as f:
        for row in csv.DictReader(f):
            dt = float(row["dt_days"])
            count = int(row["count"])
            total += count
            if dt > WINDOW_SIZE_DAYS:
                over_window += count
            if dt > WINDOW_STEP_DAYS:
                over_step += count
    if not total:
        return {}
    return {
        "n_branches": total,
        "pct_branches_over_window": 100 * over_window / total,
        "pct_branches_over_step": 100 * over_step / total,
    }


def peak_rss_gb(arm_dir):
    """Peak RSS from the /usr/bin/time -v wrapper, in GB."""
    stderr = arm_dir / "chronumental.stderr"
    if not stderr.exists():
        return None
    m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", stderr.read_text())
    return round(int(m.group(1)) / 1e6, 1) if m else None


def build_report(out_dir, arms):
    node_date_qc = load_module("node_date_qc")
    rows = []
    for arm in arms:
        arm_dir = out_dir / arm
        row = {"arm": arm, "flags": " ".join(ARMS[arm]) or "(none)"}
        row.update(parse_stdout(arm_dir / "chronumental.stdout", node_date_qc))
        row.update(parse_qc(arm_dir / "qc"))
        row.update(window_overruns(arm_dir / "qc"))
        wall = arm_dir / "wall_seconds.txt"
        row["hours"] = (
            round(float(wall.read_text()) / 3600, 2) if wall.exists() else None
        )
        row["peak_rss_gb"] = peak_rss_gb(arm_dir)
        if row.get("root_to_tip_rate") and row.get("fitted_rate"):
            row["rtt_over_fitted"] = round(
                row["root_to_tip_rate"] / row["fitted_rate"], 2
            )
        if row.get("root_to_tip_rate") and row.get("implied_subs_per_year"):
            row["implied_over_rtt"] = round(
                row["implied_subs_per_year"] / row["root_to_tip_rate"], 3
            )
        rows.append(row)

    import pandas as pd

    frame = pd.DataFrame(rows)
    csv_path = out_dir / "sweep_metrics.csv"
    frame.to_csv(csv_path, index=False)

    columns = [
        ("arm", "arm"),
        ("flags", "flags"),
        ("rtt_over_fitted", "RTT/fitted"),
        ("implied_over_rtt", "implied/RTT"),
        ("length_cor", "length_cor"),
        ("tip_residual_median_d", "med|res| d"),
        ("tip_residual_p95_d", "p95|res| d"),
        ("tip_residual_gate", "gate"),
        ("root_date", "root date"),
        ("root_date_inside_window", "root ok"),
        ("zero_mut_branch_time_pct", "0-mut time %"),
        ("pct_branches_over_window", ">90d %"),
        ("hours", "h"),
        ("peak_rss_gb", "GB"),
    ]
    lines = [
        "| " + " | ".join(label for _, label in columns) + " |",
        "|" + "|".join("---" for _ in columns) + "|",
    ]
    for row in rows:
        cells = []
        for key, _ in columns:
            value = row.get(key)
            if isinstance(value, float):
                value = f"{value:.2f}" if abs(value) < 1000 else f"{value:.0f}"
            cells.append("-" if value is None else str(value))
        lines.append("| " + " | ".join(cells) + " |")
    table = "\n".join(lines)
    (out_dir / "sweep_table.md").write_text(table + "\n")
    print(table)
    print(f"\nwrote {csv_path} and {out_dir / 'sweep_table.md'}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("mode", choices=["run", "qc", "report"])
    parser.add_argument("--tree-dir", type=Path, help="an existing output tree")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--arms",
        default=",".join(ARMS),
        help=f"comma-separated subset of {','.join(ARMS)}",
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=1000,
        help="held equal across arms; the pipeline's configured value",
    )
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument(
        "--residual-threshold",
        type=float,
        default=30.0,
        help="the pipeline's chronumental_max_median_abs_tip_residual_days",
    )
    args = parser.parse_args()

    arms = [a for a in args.arms.split(",") if a]
    unknown = [a for a in arms if a not in ARMS]
    if unknown:
        raise SystemExit(f"unknown arm(s): {unknown}; choose from {list(ARMS)}")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.mode in ("run", "qc") and not args.tree_dir:
        raise SystemExit(f"--tree-dir is required for mode {args.mode}")
    if args.mode == "run":
        run_arms(args.tree_dir, args.out_dir, arms, args.steps, args.jobs)
    elif args.mode == "qc":
        qc_arms(args.tree_dir, args.out_dir, arms, args.residual_threshold)
    else:
        build_report(args.out_dir, arms)


if __name__ == "__main__":
    main()
