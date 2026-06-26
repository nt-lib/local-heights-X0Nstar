#!/usr/bin/env bash
set -e
mkdir -p logs

shopt -s nullglob
mapfile -t magma_files < <(find computations -name '*.m' | sort)
mapfile -t sage_files  < <(find computations -name '*.sage' | sort)

if [ ${#magma_files[@]} -gt 0 ]; then
    printf '%s\n' "${magma_files[@]}" | \
        parallel --joblog logs/magma_verify_joblog.txt -j "${1:-10}" \
        'dir={//}; subdir=logs/${dir#computations/}; mkdir -p "$subdir"; magma -n {} >| "$subdir"/{/.}.txt'
fi

if [ ${#sage_files[@]} -gt 0 ]; then
    printf '%s\n' "${sage_files[@]}" | \
        parallel --joblog logs/sage_verify_joblog.txt -j "${1:-10}" \
        'dir={//}; subdir=logs/${dir#computations/}; mkdir -p "$subdir"; command time -o "$subdir"/{/.}.timing.txt -v sage {} >| "$subdir"/{/.}.txt'
fi
