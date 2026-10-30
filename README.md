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
- Optional: git 2.30 or higher, for filtering `git diff` output with `-I` (see [Checking the logs after a rerun](#checking-the-logs-after-a-rerun))

## Running computations

The code depends on several git submodules (see [Project layout](#project-layout)), so clone
the repository with
```bash
git clone --recurse-submodules https://github.com/nt-lib/local-heights-X0Nstar
cd local-heights-X0Nstar
```
In an existing clone the submodules can be initialised with
```bash
git submodule update --init
```
`make verify` initialises any submodules that are missing, and `make test` initialises
`tests/magma-unittest`.

All commands below should be run from the **repo root**:

```bash
make verify           # uses 10 parallel jobs and `magma` from the PATH by default
make verify JOBS=4 MAGMA=~/magma-2.29-10/magma
```
`make test` also accepts `MAGMA=...`.

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
The file `logs/magma_verify_joblog.txt` is the GNU parallel job log, with one line per completed job.

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

Note: some computations take many hours, see [Hardware and timings](#hardware-and-timings).

## Checking the logs after a rerun

The logs in `logs/` are committed, and a rerun overwrites them.  After the rerun has finished,
```bash
git diff
```
shows the differences between our logs and those of the rerun.  Some lines always differ: the
Magma version, date and random seed on the first line, the total time and memory usage on the
last line, and the timings that QCMod prints along the way.  For example:

```diff
--- a/logs/general/dualgraphs.txt
+++ b/logs/general/dualgraphs.txt
@@ -1,5 +1,5 @@
-Magma V2.29-10    Tue Oct 13 2026 01:21:53 on mordell  [Seed = 3228036456]
-Type ? for help.  Type <Ctrl>-D to quit.
+Magma V2.30-1     Fri Oct 30 2026 14:20:29 on mordell  [Seed = 312278796]
+Type <Ctrl>-D to quit.
@@ -57,4 +57,4 @@
-Total time: 4.490 seconds, Total memory usage: 32.09MB
+Total time: 6.530 seconds, Total memory usage: 32.09MB
```

Optionally, these lines can be hidden with the `-I<regex>` option of `git diff` (git 2.30 or
higher), which ignores changes in which all added and removed lines match one of the given
regexes.  This is only a convenience; plain `git diff` shows everything.
```bash
git diff -I '^Magma V2\.' -I '^Type ' -I '^Total time: ' -I '^Time: '
git diff --stat -I '^Magma V2\.' -I '^Type ' -I '^Total time: ' -I '^Time: '   # only list the files
```
Other differences do not mean that the reproduction failed; they are there to be reviewed.
For example, the model searches in `computations/X0_187star/model_*.m` print randomly chosen
correspondences, and Magma may list objects in a different order, such as the edges of the
dual graphs in `logs/general/dualgraphs.txt`.  Every script checks its claims with `assert`, so
a script that runs to the end without an error has verified its claims.

After inspecting the differences, `git restore logs` restores our logs.

## Hardware and timings

The logs in `logs/` were produced on a server at the University of Zagreb with an AMD EPYC 9175F
CPU with 16 cores (32 threads) and 384 GB of RAM, running Ubuntu 24.04.4 LTS, using Magma V2.30-1
and running all computations in parallel with `make verify JOBS=14`.  The server is shared, so
the timings depend on its load at the time.

There are two sources of timing information.

- GNU parallel records the wall-clock time of each job (column `JobRuntime`, in seconds) in
  `logs/magma_verify_joblog.txt`.  To list the jobs by running time:
  ```bash
  awk -F'\t' 'NR>1 {match($NF, /computations\/[^ ]*\.m/); printf "%10.1f  %s\n", $4, substr($NF, RSTART, RLENGTH)}' logs/magma_verify_joblog.txt | sort -n
  ```
  and to get the total running time when run sequentially and the wall-clock time of the
  parallel run:
  ```bash
  awk -F'\t' 'NR>1 {s += $4; if (!min || $3 < min) min = $3; if ($3 + $4 > max) max = $3 + $4} END {printf "sequential: %.1f h, parallel: %.1f h\n", s/3600, (max-min)/3600}' logs/magma_verify_joblog.txt
  ```
- Magma prints the CPU time and the memory usage of each script on the last line of its log,
  e.g. `Total time: 4.490 seconds, Total memory usage: 32.09MB`.  To list them all:
  ```bash
  grep -r "Total time" logs
  ```

## Project layout

| Path | Contents |
|------|----------|
| `src/` | Shared library code: the Magma package `local_heights.spec` (star-quotient measures, genera and bounds), attached with `AttachSpec("local_heights.spec")`, and patched copies of QCMod files (Mordell-Weil sieve, misc utilities), loaded with `load` |
| `computations/X0_<N>star/` | Scripts for the curve X_0(N)^* |
| `computations/general/` | Cross-cutting verifications: existence of clean correspondences, fixed-point formulas, dual graphs, etc. |
| `logs/` | Output committed as a record of past computations |
| `data/` | Precomputed data files (e.g. class numbers) |
| `tests/` | Unit tests for `src/`, run with `make test` |
| `QCMod/` | Submodule: quadratic Chabauty engine |
| `ModularCurvesX0plusG4-6/` | Submodule: model-finding code for X_0(N)^* |
| `ShimuraCurveALQuotients/` | Submodule: rational CM points on Atkin-Lehner quotients (Assaf-Hashimoto) |
| `tests/magma-unittest/` | Submodule: unit test framework, only needed for `make test` |

Within each `computations/X0_<N>star/` directory the naming convention is:

| Prefix | Meaning |
|--------|---------|
| `qc_p<prime>.m` | Full quadratic Chabauty run at the given prime |
| `qc_p<prime>_mws.m` | QC run followed by a Mordell-Weil sieve step |
| `model_p<prime>.m` | Model search via `model_equation_finder` at the given prime |
| `model_search_*.m` | Ad-hoc / random model search |

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
