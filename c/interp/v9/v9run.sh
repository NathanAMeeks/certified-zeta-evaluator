#!/bin/bash
# v9 driver (combination of v7 speed + v8 self-correction, directed by Nathan A. Meeks; v8 gated self-correcting loop credited to him).
# usage: v9run.sh T P B   (window args as v7; AHCACHE auto-extended by the v8 delta extender in fftzi_v9)
# Attempt 1 = v7 settings exactly (delta 0.24, caller's P). Only if it is not PROVEN (Lemma CP) does the ladder escalate:
#   same delta with P doubled until P*delta >= 1.2*(B+140), then delta 0.2, 0.16, 0.12 (P grown the same way).
# Every attempt is an independent full certification; only a PROVEN result is accepted. Ladder is search-only.
D=${V9DIR:-$(cd "$(dirname "$0")/.." && pwd)}; T=$1; P=$2; B=$3; export AHCACHE=${AHCACHE:-$D/ahcache.bin}
need=$(( (12*(B+140)+9)/10 )); spent=0; o=""
run(){ o=$(taskset -c ${CPU:-7} $D/fftzi_v9 $T $1 $2 80 $B 60 /dev/null dbl 2>/dev/null)
  s=$(echo "$o" | grep -o '"sec_total":[0-9.]*' | cut -d: -f2); spent=$(awk -v a=$spent -v b=${s:-0} 'BEGIN{print a+b}')
  if echo "$o" | grep -q '"status":"PROVEN"'; then echo "$o" | sed "s/}\$/,\"v9_delta\":$1,\"v9_P\":$2,\"sec_all_attempts\":$spent}/"; exit 0; fi
  echo "v9run: T=$T delta=$1 P=$2 not proven, escalate" >&2; }
grow(){ local Pn=$1 dl=$2; while awk -v p=$Pn -v d=$dl -v n=$need 'BEGIN{exit !(p*d<n)}'; do Pn=$((Pn*2)); done; echo $Pn; }
run 0.24 $P
Pg=$(grow $P 0.24); [ "$Pg" != "$P" ] && run 0.24 $Pg
for dl in 0.2 0.16 0.12; do run $dl $(grow $P $dl); done
echo "$o"; exit 1
