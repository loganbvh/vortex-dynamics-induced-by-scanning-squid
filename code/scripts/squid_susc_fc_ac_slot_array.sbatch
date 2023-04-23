#!/bin/bash
#SBATCH --job-name=squid_susc_fc_ac_slot_array
#SBATCH --array=0-100
#SBATCH --output=slurm-%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8GB
#SBATCH --time=24:00:00

set -x

# Activate conda env
source $GROUP_HOME/miniconda3/bin/activate
conda activate tdgl

export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}

outdir=$GROUP_SCRATCH/py-tdgl/results/$SLURM_ARRAY_JOB_ID
mkdir -p $outdir

pyscript=$HOME/nb-near-tc/scripts/squid_susc_fc_ac.py

# Copy the python script and this shell script to the results directory
cp -u $pyscript $outdir/
cp -u $0 $outdir/

python $pyscript \
    --directory=$outdir/$SLURM_ARRAY_TASK_ID \
    --I_fc 2 \
    --slot-size 2 50 \
    --slot-radius=0.1 \
    --slot-top-center-x -6 6 $SLURM_ARRAY_TASK_COUNT \
    --slot-top-center-y 0 0 $SLURM_ARRAY_TASK_COUNT \
    --index=$SLURM_ARRAY_TASK_ID \
    --squid-type=hypres-small \
    --squid-position 0 0 0.5 \
    --squid-points=4000 \
    --squid-smooth=100 \
    --squid-iterations=5 \
    --film-radius=20 \
    --film-points=4000 \
    --film-smooth=100 \
    --max-edge-length=0.5 \
    --cycles=1.5 \
    --points-per-cycle=51 \
    --d=0.2 \
    --lam=1.2 \
    --xi=0.9 \
    --gamma=1 \
    --solve-time 20e3 2e3 \
    --dt-init=1e-6 \
    --save-every=500 \
    --seed-solutions

# Move the stdout log to the results directory
mv "slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" $outdir
