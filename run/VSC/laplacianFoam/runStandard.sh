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


# OpenFOAM standard

#Use the controlDict with the correct OptimizationSwitch
mv system/controlDict system/controlDictComm
mv system/controlDictNormal system/controlDict
mpirun -np ${n} laplacianFoam -parallel > log.foamRun
mv system/controlDict system/controlDictNormal
mv system/controlDictComm system/controlDict
