# Contributing to SMVac-Research

## Development setup

- C++17 compiler (GCC, Clang, or MSVC), CMake >= 3.10.
- Python 3 with `numpy`, `matplotlib`, `pandas` (only for the scan driver
  and figure scripts).

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

## Repository rules that keep the physics honest

- **Any intentional change to the physics or to production numerical
  settings requires regenerating the frozen references**: run
  `./build/write_reference_tables`, inspect the diff, and record the old
  and new benchmark values in `docs/audit-notes.md` and this changelog.
  The regression tests fail by design after such a change; that is the
  signal to regenerate consciously, never to bypass.
- `tests/physics/` contains the validations that back the scientific
  claims (analytic limits, convergence, literature comparison); extend
  them when adding capabilities, and update `docs/validation.md` with the
  measured numbers.
- Keep the documentation honest: claims must be backed by a test or a
  committed dataset, and limitations belong in `docs/limitations.md`.

## Scan and figure workflow

```bash
python scripts/run_phase_diagram.py --mode numerical \
    --mt-min 155 --mt-max 185 --mh-min 112 --mh-max 138 --step 0.25 \
    --output data/numerical_sm_region.csv
python scripts/make_figures.py
```

`results/` is scratch space (gitignored); committed datasets live in
`data/` with a README describing exactly how each was produced.

## Style

- The physics modules (`src/`) favor transparency over abstraction: a
  reader should be able to map equations in `docs/theory.md` to lines of
  code one-to-one (`docs/numerical-method.md` provides the mapping table).
- Comments explain conventions and constraints (units, loop-counting
  factors, the `t = ln(mu^2)` evolution variable), not the obvious.
