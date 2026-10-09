#!/usr/bin/env bash
# Optional: (re)writes data/classnumbers_upto_<N_max>.dat, the class number tables read by
# src/genera_and_bounds.m, from the LMFDB tables of class groups of imaginary quadratic fields
# (https://www.lmfdb.org/NumberField/QuadraticImaginaryClassGroups, k = 0).  The .dat file is
# part of the repository, so this is only needed to check or rebuild it; verify_all.sh does
# not run it.  The four LMFDB files take about 1.1 GB uncompressed.
set -e
cd "$(dirname "$0")"
mkdir -p data

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

magma -n initialize_class_numbers.m < /dev/null
