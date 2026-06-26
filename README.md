# Vanishing Local Heights on $X_0(N)^*$

This repository is about studying rational points on $X_0(N)^*$, the quotient
of the modular curve $X_0(N)$ by all Atkin-Lehner involutions.  We develop a
general method for constructing correspondences for quadratic Chabauty whose
local height contributions at the primes of bad reduction are trivial.  For
executing quadratic Chabauty we rely on the
[QCMod](https://github.com/steffenmueller/QCMod) package.

As an application of our methods we determine the rational points on
$X_0(187)^{\ast}$, $X_0(247)^{\ast}$, and $X_0(319)^{\ast}$.

## Requirements

- Magma v2.28 or higher
- GNU parallel (for `make verify`)

## Submodule dependencies

| Submodule | Role |
|-----------|------|
| `QCMod/` | Quadratic Chabauty engine |
| `ModularCurvesX0plusG4-6/` | Model-finding code for X_0(N)^* |

Initialise before running any computations:

```bash
git submodule update --init
```

## Project layout

| Path | Contents |
|------|----------|
| `src/` | Shared library code: star-quotient measures, rank computations, Mordell-Weil sieve, misc utilities |
| `computations/X0_<N>star/` | Scripts for the curve X_0(N)^* |
| `computations/general/` | Cross-cutting verifications: genus tables, fixed-point formulas, dual graphs, etc. |
| `logs/` | Output committed as a record of past computations |
| `notebooks/` | Exploratory Jupyter notebooks |
| `data/` | Precomputed data files (e.g. class numbers) |

Within each `computations/X0_<N>star/` directory the naming convention is:

| Prefix | Meaning |
|--------|---------|
| `qc_p<prime>.m` | Full quadratic Chabauty run at the given prime |
| `qc_p<prime>_mws.m` | QC run followed by a Mordell-Weil sieve step |
| `model_p<prime>.m` | Model search via `model_equation_finder` at the given prime |
| `model_search_*.m` | Ad-hoc / random model search |
| `gonality.m` | Gonality / gonal map computation |
| `reductions.m` | Reduction checks |

## Running computations

All commands below should be run from the **repo root**:

```bash
make verify           # uses 10 parallel jobs by default
```

Alternatively one can run the `verify_all.sh` script manually using
```bash
./verify_all.sh 4     # the 4 means use 4 parallel jobs
```

Output from each script is captured to `logs/<curve>/<script>.txt`, mirroring the `computations/` layout.

Note: QC computations (especially for X_0(319)^*) can take many hours.

## Running on a remote machine

```bash
make verify_remote ssh="user@hostname"  # sync and run all computations
make copy_logs     ssh="user@hostname"  # fetch logs back afterwards
```

## Troubleshooting

**`parallel: command not found`**

GNU parallel is not installed.  Install it with:

```bash
# macOS
brew install parallel

# Debian / Ubuntu
sudo apt install parallel

# RHEL / Fedora
sudo dnf install parallel
```
