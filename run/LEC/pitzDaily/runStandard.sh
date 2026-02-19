#!/bin/bash
#SBATCH --job-name=pitz_2                 # Job name
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


# OpenFOAM standard

#Use the controlDict with the correct OptimizationSwitch
mv system/controlDict system/controlDictComm
mv system/controlDictNormal system/controlDict
mpirun --mca btl_openib_allow_ib true -np ${n} foamRun -parallel > log.foamRun
mv system/controlDict system/controlDictNormal
mv system/controlDictComm system/controlDict
