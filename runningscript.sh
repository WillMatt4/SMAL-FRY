#!/bin/sh
##up to 12h:
#SBATCH --partition=shared-cpu
##12h - 4days:
##SBATCH --partition=public-cpu
##up to 15min:
##SBATCH --partition=debug-cpu
##max. of 2 cores for 14 days:
##SBATCH --partition=public-longrun-cpu

##SBATCH --exclude=cpu[116-119]

#SBATCH --ntasks=5
#SBATCH --cpus-per-task=1  
##SBATCH --time=24:00:00

##SBATCH --mail-type=ALL

#SBATCH -J deltazeta_vanilla##test0##rdn0_1000
#SBATCH --output=slurm-%x-%J.out
e

echo "Starting at `date`"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Current working directory is `pwd`"
echo ""
echo "***** LAUNCHING *****"
echo `date '+%F %H:%M:%S'`
echo ""

# load Anaconda and OpenMPI
#module load GCCcore/6.4.0
#module load glibc/2.17
#module load GCCcore/8.3.0
#module load glibc/2.30
module load Anaconda3
#module load GCC/11.3.0
#module load OpenMPI/4.1.4
module load foss

echo "Loaded Anaconda3 and foss"
echo ""

# # #OpenMP settings:
# export OMP_NUM_THREADS=2
# export OMP_PLACES=threads
# export OMP_PROC_BIND=spread



#launch simulations
srun python3 SMAL-FRY.py

echo ""
echo "***** DONE *****"
echo `date '+%F %H:%M:%S'`
echo ""
