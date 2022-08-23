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

#SBATCH --ntasks=8 
## 5 params + fiducial + 2 bias 
##SBATCH --cpus-per-task=3  
#SBATCH --time=12:00:00

#SBATCH --mail-type=ALL

#SBATCH -J deltazeta##test0##rdn0_1000
#SBATCH --output=slurm-%x-%J.out


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
#module load Anaconda3
##module load GCC/10.2.0
##module load OpenMPI/4.0.5
#module load foss
#conda activate myenv
##module load mpi4py/3.0.3-timed-pingpong
#module load numpy/1.18.5-Python-3.6.6
#module load scipy/1.4.1-Python-3.7.4

#echo "Loaded Anaconda3 and foss"
echo ""

# # #OpenMP settings:
# export OMP_NUM_THREADS=2
# export OMP_PLACES=threads
# export OMP_PROC_BIND=spread



#launch simulations
srun python3 SMAL-FRY.py
wait
echo ""
echo "***** DONE *****"
echo `date '+%F %H:%M:%S'`
echo ""
