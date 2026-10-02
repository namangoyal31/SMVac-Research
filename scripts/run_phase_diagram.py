#!/usr/bin/env python
"""Parallel parameter-scan driver for the SM vacuum stability classifiers.

Enumerates the (Mh, Mt) grid exactly like apps/generate_phase_diagram.cpp
(Mt-major), splits it into index-range chunks, runs the compiled binary on
each chunk in parallel, and aggregates the chunk CSVs into a single file.

Example:
    python scripts/run_phase_diagram.py --mode numerical \
        --mt-min 160 --mt-max 185 --mh-min 110 --mh-max 140 --step 0.5

The binary is expected at build/generate_phase_diagram (override with
--binary). Output columns and status codes are documented in
docs/numerical-method.md.
"""

import argparse
import concurrent.futures
import glob
import os
import subprocess
import sys

import pandas as pd


def build_grid(args):
    """Mt-major enumeration matching the C++ app."""
    n_mt = round((args.mt_max - args.mt_min) / args.step) + 1
    n_mh = round((args.mh_max - args.mh_min) / args.step) + 1
    mts = [args.mt_min + args.step * i for i in range(n_mt)]
    mhs = [args.mh_min + args.step * j for j in range(n_mh)]
    return [(mt, mh) for mt in mts for mh in mhs]


def run_chunk(binary, args, start, end, chunk_dir, retries=2):
    produced = os.path.join(chunk_dir, f"{args.mode}_chunk_{start}.csv")
    cmd = [
        binary,
        f"--{args.mode}",
        "--mt-min", str(args.mt_min),
        "--mt-max", str(args.mt_max),
        "--mh-min", str(args.mh_min),
        "--mh-max", str(args.mh_max),
        "--step", str(args.step),
        "--start", str(start),
        "--end", str(end),
        "--output-dir", chunk_dir,
        "--precision", str(args.precision),
        "--c6", str(args.c6),
    ]
    # Retry: transient file locks (e.g. antivirus scans on Windows) can make
    # a chunk fail spuriously; a failed chunk is rerun from scratch.
    last_err = None
    for attempt in range(retries + 1):
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0 and os.path.isfile(produced):
            return produced
        last_err = (result.returncode, result.stdout, result.stderr)
    sys.stderr.write(last_err[1])
    sys.stderr.write(last_err[2])
    raise RuntimeError(f"chunk [{start}, {end}) failed after "
                       f"{retries + 1} attempts (code {last_err[0]})")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mode", choices=["analytical", "numerical"], default="numerical")
    parser.add_argument("--mt-min", type=float, default=160.0)
    parser.add_argument("--mt-max", type=float, default=185.0)
    parser.add_argument("--mh-min", type=float, default=110.0)
    parser.add_argument("--mh-max", type=float, default=140.0)
    parser.add_argument("--step", type=float, default=0.5)
    parser.add_argument("--processes", type=int, default=os.cpu_count() or 4)
    parser.add_argument("--binary",
                        default=os.path.join("build", "generate_phase_diagram" + (".exe" if os.name == "nt" else "")),
                        help="path to the compiled generate_phase_diagram binary")
    parser.add_argument("--output", default=None,
                        help="aggregated output CSV (default: results/<mode>_data.csv)")
    parser.add_argument("--precision", type=int, default=12,
                        help="significant digits in the CSV output")
    parser.add_argument("--c6", type=float, default=1.0,
                        help="coefficient of the Planck-suppressed phi^6 operator")
    parser.add_argument("--results-dir", default="results")
    parser.add_argument("--keep-chunks", action="store_true")
    args = parser.parse_args()

    output = args.output or os.path.join(args.results_dir, f"{args.mode}_data.csv")

    grid = build_grid(args)
    total = len(grid)
    print(f"Grid: {total} points (Mt in [{args.mt_min}, {args.mt_max}], "
          f"Mh in [{args.mh_min}, {args.mh_max}], step {args.step})")

    if not os.path.isfile(args.binary):
        sys.exit(f"Binary not found at '{args.binary}'. Build first, e.g.:\n"
                 "  cmake -S . -B build && cmake --build build")

    os.makedirs(args.results_dir, exist_ok=True)
    chunk_dir = os.path.join(args.results_dir, "chunks")
    os.makedirs(chunk_dir, exist_ok=True)

    n_chunks = max(1, args.processes)
    chunk_size = (total + n_chunks - 1) // n_chunks
    ranges = [(s, min(s + chunk_size, total)) for s in range(0, total, chunk_size)]

    print(f"Running {len(ranges)} chunks on {args.processes} workers...")
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.processes) as pool:
        futures = [pool.submit(run_chunk, args.binary, args, s, e, chunk_dir)
                   for (s, e) in ranges]
        chunk_files = [f.result() for f in concurrent.futures.as_completed(futures)]

    df = pd.concat((pd.read_csv(f) for f in sorted(chunk_files)), ignore_index=True)
    df.to_csv(output, index=False)
    print(f"Aggregated {len(df)} rows -> {output}")

    if not args.keep_chunks:
        for f in glob.glob(os.path.join(chunk_dir, "*.csv")):
            os.remove(f)


if __name__ == "__main__":
    main()
