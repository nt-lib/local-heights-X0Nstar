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
| `ShimuraCurveALQuotients/` | Rational CM points on Atkin-Lehner quotients (Assaf-Hashimoto) |

`make verify` initialises any submodules that are missing, and `make test` initialises
`tests/magma-unittest`. To initialise them by hand:

```bash
git submodule update --init
```

## Project layout

| Path | Contents |
|------|----------|
| `src/` | Shared library code: the Magma package `local_heights.spec` (star-quotient measures, genera and bounds), attached with `AttachSpec("local_heights.spec")`, and patched copies of QCMod files (Mordell-Weil sieve, misc utilities), loaded with `load` |
| `computations/X0_<N>star/` | Scripts for the curve X_0(N)^* |
| `computations/general/` | Cross-cutting verifications: existence of clean correspondences, fixed-point formulas, dual graphs, etc. |
| `logs/` | Output committed as a record of past computations |
| `data/` | Precomputed data files (e.g. class numbers) |

Within each `computations/X0_<N>star/` directory the naming convention is:

| Prefix | Meaning |
|--------|---------|
| `qc_p<prime>.m` | Full quadratic Chabauty run at the given prime |
| `qc_p<prime>_mws.m` | QC run followed by a Mordell-Weil sieve step |
| `model_p<prime>.m` | Model search via `model_equation_finder` at the given prime |
| `model_search_*.m` | Ad-hoc / random model search |

## Running computations

All commands below should be run from the **repo root**:

```bash
make verify           # uses 10 parallel jobs by default
```

Alternatively one can run the `verify_all.sh` script manually using
```bash
./verify_all.sh 4     # the 4 means use 4 parallel jobs
```
An optional second argument selects the Magma executable (default: `magma` from the `PATH`):
```bash
./verify_all.sh 4 ~/magma-2.29-10/magma
```
Magma refuses to load package files (`*.sig`) compiled by a newer Magma version; if you switch to
an older Magma, run `make clean` first (see Troubleshooting below).

Output from each script is captured to `logs/<curve>/<script>.txt`, mirroring the `computations/` layout.

`computations/general/existence_of_relations.m` reads class numbers from
`data/classnumbers_upto_2589510.dat`, which is part of this repository.  This file was written
from the LMFDB data files (about 1.1 GB uncompressed, not part of this repository) by
```bash
./initialize_class_numbers.sh
```
which you can run to check or rebuild it (`verify_all.sh` does not run it).  It downloads the
LMFDB files to `data/` if they are missing (this needs `curl` and `gunzip`); alternatively,
download the four files ($|d| \equiv 3 \pmod 8$, $7 \pmod 8$, $4 \pmod{16}$,
$8 \pmod{16}$, with $k = 0$ for all four) from
<https://www.lmfdb.org/NumberField/QuadraticImaginaryClassGroups> and uncompress them into
`data/` as `cl3mod8.0`, `cl7mod8.0`, `cl4mod16.0`, `cl8mod16.0`.

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

**`Package ... is compiled for newer version than current executable`, or `Undefined reference`
errors inside a package (e.g. in `ShimuraCurveALQuotients/`)**

Magma caches compiled packages as `*.sig` files next to the sources.  These can be stale, e.g.
after switching to an older Magma version or after checking out a different version of a
submodule.  Delete all of them, including those in the submodules, with
```bash
make clean
```
and run the computations again; Magma recompiles the packages as needed.  Do not do this while
computations are running in the same checkout.

## License

This repository is licensed under the GNU General Public License v3.0; see [LICENSE](LICENSE).
The submodules (QCMod, ModularCurvesX0plusG4-6 and magma-unittest) are separate projects with their own terms.
