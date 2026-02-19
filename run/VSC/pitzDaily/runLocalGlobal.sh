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


# With local and global communication optimization

mpirun -np ${n} optimizeCommPar -parallel  > log.optimizeCommPar
mpirun -rf ./constant/rankFile -np ${n} foamRun -parallel > log.foamRun
