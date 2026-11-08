#!/bin/bash
R=/workspace/rh_extend; S=$R/scratch/t106c; NAME=$1; shift; B=/workspace/audit_scratch/t106c/$NAME
rm -rf $B; mkdir -p $B/brackets $B/logs $B/scratch; cp $R/ahcache.bin $B/; cp $R/scratch/t102/bench_1core/ck0.json $B/checkpoint.json
L0=$(cut -d' ' -f1 /proc/loadavg)
cd $B && env WORKER_CPUS=1,2,3,4,5,6,7 GATE_CPUS=0,1,2,3,4,5,6,7 GATE_PROCS=14 GATE_K=8 XZCPUS=0-7 XZPAR=6 CPU=0 RHX_HOME=$B BATCH_WINDOWS=400 T102=1 GTURING_P=1 HYB_VC=0 HYB_TMESH=1 HYB_REUSE=0 CK_EVERY=10 GRS_TIER1=1 GRS_BIN_DIR=$R/spotcheck/grs "$@" nice -n 19 python3 $S/run.py python3 $S/driver_t106.py
echo "$NAME load=$L0-$(cut -d' ' -f1 /proc/loadavg) $(cat time.txt)" >> /workspace/audit_scratch/t106c/res.txt
