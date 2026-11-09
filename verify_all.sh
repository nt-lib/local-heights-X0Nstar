#!/usr/bin/env bash
# Usage: ./verify_all.sh [jobs] [magma]
#   jobs   number of parallel jobs (default 10)
#   magma  Magma executable to use (default: magma from the PATH)
set -e

export MAGMA="${2:-magma}"
if ! command -v "$MAGMA" > /dev/null; then
    echo "Error: Magma executable '$MAGMA' not found" >&2
    exit 1
fi

# The computations load code from the submodules (QCMod, ModularCurvesX0plusG4-6, ...).
# Without them every such script fails, but Magma still exits with status 0.
missing=$(git submodule status | awk '/^-/ {print $2}')
if [ -n "$missing" ]; then
    echo "Error: the following submodules are not checked out:" >&2
    echo "$missing" | sed 's/^/    /' >&2
    echo "Run the following command from the repository root and try again:" >&2
    echo "    git submodule update --init" >&2
    exit 1
fi

mkdir -p logs

shopt -s nullglob
mapfile -t magma_files < <(find computations -name '*.m' | sort)
mapfile -t sage_files  < <(find computations -name '*.sage' | sort)

if [ ${#magma_files[@]} -gt 0 ]; then
    printf '%s\n' "${magma_files[@]}" | \
        parallel --joblog logs/magma_verify_joblog.txt -j "${1:-10}" \
        'dir={//}; subdir=logs/${dir#computations/}; mkdir -p "$subdir"; "$MAGMA" -n {} >| "$subdir"/{/.}.txt'
fi

if [ ${#sage_files[@]} -gt 0 ]; then
    printf '%s\n' "${sage_files[@]}" | \
        parallel --joblog logs/sage_verify_joblog.txt -j "${1:-10}" \
        'dir={//}; subdir=logs/${dir#computations/}; mkdir -p "$subdir"; command time -o "$subdir"/{/.}.timing.txt -v sage {} >| "$subdir"/{/.}.txt'
fi
