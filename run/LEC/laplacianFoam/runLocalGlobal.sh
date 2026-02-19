#!/bin/bash
#SBATCH --job-name=laplacianFoam_2        # Job name
#SBATCH --nodes=2                         # Number of nodes
#SBATCH --ntasks-per-node=48              # Number of MPI tasks per node
#SBATCH --cpus-per-task=1                 # Number of CPU cores per task
#SBATCH --time=01:00:00                   # Walltime
#SBATCH --partition=hpc2                  # Queue/Partition name
#SBATCH --account=sieberer                # Account name
#SBATCH --output=openfoam_job_%j.out      # Output file
#SBATCH --error=openfoam_job_%j.err       # Error file
#SBATCH --wckey=....

source /home/sieberer/OpenFOAM/OpenFOAM-12/etc/bashrc

nn=2
nc=48
n=$(( nn * nc ))



# With local and global communication optimization

mpirun --mca btl_openib_allow_ib true -np ${n} optimizeCommPar -parallel  > log.optimizeCommPar

mpirun --mca btl_openib_allow_ib true -rf ./constant/rankFile -np ${n} laplacianFoam -parallel > log.foamRun
