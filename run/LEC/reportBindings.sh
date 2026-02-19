#!/bin/bash
#SBATCH --job-name=reportBindings         # Job name
#SBATCH --nodes=1                         # Number of nodes
#SBATCH --ntasks-per-node=48              # Number of MPI tasks per node
#SBATCH --cpus-per-task=1                 # Number of CPU cores per task
#SBATCH --time=00:10:00                   # Walltime
#SBATCH --partition=hpc2                  # Queue/Partition name
#SBATCH --account=sieberer                # Account name
#SBATCH --output=openfoam_job_%j.out      # Output file
#SBATCH --error=openfoam_job_%j.err       # Error file
#SBATCH --wckey=....

source /home/sieberer/OpenFOAM/OpenFOAM-12/etc/bashrc

mpirun --report-bindings -bind-to core -np 48 /bin/true  > reportBindings.txt 2>&1

