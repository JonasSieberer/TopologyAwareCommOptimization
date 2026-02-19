#!/bin/bash
#SBATCH --job-name=OF_prep                # Job name
#SBATCH --nodes=1                         # Number of nodes
#SBATCH --ntasks-per-node=1               # Number of MPI tasks per node
#SBATCH --cpus-per-task=1                 # Number of CPU cores per task
#SBATCH --time=1:00:00                    # Walltime
#SBATCH --partition=hpc2                  # Queue/Partition name
#SBATCH --account=sieberer                # Account name
#SBATCH --output=openfoam_job_%j.out      # Output file
#SBATCH --error=openfoam_job_%j.err       # Error file
#SBATCH --wckey=XT11.2

source /home/sieberer/OpenFOAM/OpenFOAM-12/etc/bashrc

blockMesh > log.blockMesh

decomposePar -force -cellProc > log.decomposePar


