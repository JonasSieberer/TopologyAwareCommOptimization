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


# Only global communication optimization

# Copy the reportBindings.txt file into the case directory
cp ../reportBindings.txt ./reportBindings.txt

# Use the decompositionDict with the version 2 for optimizeCommPar
mv system/decomposeParDict system/decomposeParDictCorr
mv system/decomposeParDict2 system/decomposeParDict
mpirun --mca btl_openib_allow_ib true -np ${n} optimizeCommPar -parallel  > log.optimizeCommPar
mv system/decomposeParDict system/decomposeParDict2
mv system/decomposeParDictCorr system/decomposeParDict

mpirun --mca btl_openib_allow_ib true -bind-to core -np ${n} foamRun -parallel > log.foamRun
