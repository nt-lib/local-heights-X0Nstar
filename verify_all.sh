#!/usr/bin/env bash
set -e
mkdir -p logs data

# Class group data from the LMFDB (k = 0), read by computations/general/existence_of_relations.m.
# A valid data/cl*.0.gz (e.g. downloaded in a browser) is used instead of downloading.
for f in cl3mod8 cl7mod8 cl4mod16 cl8mod16; do
    if [ ! -f "data/$f.0" ]; then
        if ! gzip -t "data/$f.0.gz" 2>/dev/null; then
            # Keep a partial download only if it starts like a gzip file (not an HTML error page).
            if [ "$(head -c 2 "data/$f.0.gz" 2>/dev/null | od -An -tx1 | tr -d ' \n')" != 1f8b ]; then
                rm -f "data/$f.0.gz"
            fi
            # Up to 7 tries, 30 s apart; -C - continues where the previous try stopped.
            for try in 1 2 3 4 5 6 7; do
                if curl -fL -C - -o "data/$f.0.gz" "https://www.lmfdb.org/NumberField/QuadraticImaginaryClassGroups?filenamebase=$f&k=0&Fetch=fetch"; then
                    break
                fi
                if [ "$try" -lt 7 ]; then
                    echo "Download of $f.0.gz failed (try $try of 7); continuing in 30 s." >&2
                    sleep 30
                fi
            done
        fi
        gunzip "data/$f.0.gz" || {
            echo "Could not unpack data/$f.0.gz. If the LMFDB blocked the download, get $f.0.gz (k = 0) in a browser from https://www.lmfdb.org/NumberField/QuadraticImaginaryClassGroups and put it in data/." >&2
            exit 1
        }
    fi
done

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
