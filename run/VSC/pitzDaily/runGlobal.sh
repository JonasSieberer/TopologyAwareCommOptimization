#!/bin/bash
#SBATCH -A p.....
#SBATCH -J OF_run
#SBATCH -t 00:20:00
#SBATCH --partition=zen3_0512
#SBATCH --qos=zen3_0512
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128

source $HOME/OpenFOAM/OpenFOAM-12/etc/bashrc

nn=2
nc=128
n=$(( nn * nc ))

# Only global communication optimization

# Copy the reportBindings.txt file into the case directory
cp ../reportBindings.txt ./reportBindings.txt

# Use the decompositionDict with the version 2 for optimizeCommPar
mv system/decomposeParDict system/decomposeParDictCorr
mv system/decomposeParDict2 system/decomposeParDict
mpirun -np ${n} optimizeCommPar -parallel  > log.optimizeCommPar
mv system/decomposeParDict system/decomposeParDict2
mv system/decomposeParDictCorr system/decomposeParDict

mpirun -bind-to core -np ${n} foamRun -parallel > log.foamRun
