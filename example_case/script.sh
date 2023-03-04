#!/bin/bash
#SBATCH --partition=workq
#SBATCH --nodes=4
#SBATCH --time=24:00:00
#SBATCH --job-name="Red_p01_c0.6"
#SBATCH -A kXXX
#SBATCH --err=std.err
#SBATCH --output=std.out

echo "start:"
echo "The job" ${SLURM_JOB_ID} "is running on" ${SLURM_JOB_NODELIST}
date
#
# generating mesh
# time srun -N 1 blockMesh
# partioning mesh
time srun -N 1 decomposePar -fileHandler collated

time srun --ntasks=128 --hint=nomultithread --ntasks-per-node=32 --ntasks-per-socket=16 your_solver_name -parallel -fileHandler collated

#time srun -N 1 reconstructPar -latestTime
#rm -rf proces*
# cleaning up
#tar -cf .tar *

echo "finish:"
date
