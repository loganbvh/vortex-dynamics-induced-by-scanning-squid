#!/bin/bash

outdir="../../data/simulation/squid-single-vortex"

if [[ -e $outdir || -L $outdir ]] ; then
    i=0
    while [[ -e $outdir-$i || -L $outdir-$i ]] ; do
        let i++
    done
    outdir=$outdir-$i
fi

mkdir -p $outdir


python squid_susc_fc_ac.py \
    --directory=$outdir \
    --I_fc 2.5 2.5 1 \
    --index=0 \
    --squid-type=hypres-small \
    --squid-angle=0 \
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
    --solve-time 500 500 \
    --dt-init=1e-6 \
    --save-every=100 \
    --seed-solutions
