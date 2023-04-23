#!/bin/bash

outdir="../../data/simulation/squid-slot-1"

slot_top_center_x=$1
slot_top_center_y=$2

outdir=$outdir/"xy_${slot_top_center_x}_${slot_top_center_y}"
mkdir -p $outdir


python squid_susc_fc_ac.py \
    --directory=$outdir \
    --I_fc 1.0 \
    --slot-size 2 50 \
    --slot-radius=0.1 \
    --slot-top-center-x $slot_top_center_x $slot_top_center_x 1 \
    --slot-top-center-y $slot_top_center_y $slot_top_center_y 1 \
    --index=0 \
    --squid-type=hypres-small \
    --squid-angle=180 \
    --squid-position 0 0 0.5 \
    --squid-points=4000 \
    --squid-smooth=100 \
    --squid-iterations=5 \
    --film-radius=15 \
    --film-shape=box \
    --film-points=4000 \
    --film-smooth=100 \
    --cycles=0 \
    --points-per-cycle=8 \
    --d=0.2 \
    --lam=1.35 \
    --xi=0.9 \
    --gamma=1 \
    --solve-time 5e3 2e3 \
    --dt-init=1e-6 \
    --save-every=500 \
    --seed-solutions
