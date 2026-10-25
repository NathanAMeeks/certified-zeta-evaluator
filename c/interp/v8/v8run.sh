#!/bin/bash
# v8 self-correcting driver (Nathan A. Meeks' gated self-correction). usage: v8run.sh T P B  (window args as v7; AHCACHE auto-extended)
# Ladder is search-only: each attempt is an independent full certification; only a PROVEN (Lemma CP) result is accepted.
D=${V8DIR:-$(cd "$(dirname "$0")/.." && pwd)}; T=$1; P=$2; B=$3; export AHCACHE=${AHCACHE:-$D/ahcache.bin}
for dl in 0.24 0.2 0.16 0.12; do Pn=$P; while python3 -c "import sys;sys.exit(0 if $Pn*$dl<1.2*($B+140) else 1)"; do Pn=$((Pn*2)); done
  o=$(taskset -c ${CPU:-7} $D/fftzi_v8 $T $dl $Pn 80 $B 60 /dev/null dbl 2>/dev/null)
  if echo "$o" | grep -q '"status":"PROVEN"'; then echo "$o" | sed "s/}\$/,\"v8_delta\":$dl}/"; exit 0; fi
  echo "v8run: T=$T delta=$dl not proven, retry" >&2; done; echo "$o"; exit 1
