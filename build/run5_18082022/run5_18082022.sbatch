#!/bin/bash
# Job name:
#SBATCH --job-name=run5_18082022
#
# Partition:
#SBATCH --partition=CMT
#
# Processors:
#SBATCH --ntasks=32
#SBATCH --nodes=2
#
# Wall clock limit:
#SBATCH --time=24:00:00
#
# Email at end:
#SBATCH --mail-user=sxd190113@utdallas.edu
#SBATCH --mail-type=END
#
# Specify output files:
#SBATCH -o run5_18082022.o%j
#------------------------------------------------------

echo "Starting job run5_18082022"

module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=run5_18082022.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job run5_18082022 complete"
