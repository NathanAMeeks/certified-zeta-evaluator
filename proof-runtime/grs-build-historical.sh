#!/bin/bash
# [NUM] t106 G-RS batch build (candidate, unrefereed; Nathan A. Meeks & Grok). Source = scratch/grs_stage/src (Audit fixes 1-3 applied there, NOT yet re-audited).
# Pinned ISA per Audit item 4; -ffp-contract=off -fno-fast-math always.
set -e; cd "$(dirname "$0")"; F=/workspace/rh-research/tools/flintc/root/usr; L=$F/lib/x86_64-linux-gnu
CC="nice -n 19 gcc -fopenmp -I$F/include -I$F/include/x86_64-linux-gnu -Isrc"; SAFE="-ffp-contract=off -fno-fast-math"; LD="-L$L -lflint -lm -Wl,-rpath,$L"
$CC -O3 -march=x86-64-v4 -mprefer-vector-width=256 $SAFE -o bin/rs_fast src/rs_fast.c $LD
$CC -O3 -march=x86-64-v4 -mprefer-vector-width=512 $SAFE -o bin/rs_fast512 src/rs_fast.c $LD
$CC -O2 -march=x86-64-v3 $SAFE -DKTERMS=22 -DYTARGET=1.5 -DNDIRECT=1000 -o bin/rs_tb_fix src/rs_tb_fix.c $LD
sha256sum bin/rs_fast bin/rs_fast512 bin/rs_tb_fix | tee bin/SHA256
