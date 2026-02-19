#!/bin/bash
#SBATCH -A p.....
#SBATCH -J
#SBATCH -t 00:20:00
#SBATCH --partition=zen3_0512
#SBATCH --qos=zen3_0512
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128

source $HOME/OpenFOAM/OpenFOAM-12/etc/bashrc

mpirun --report-bindings -bind-to core -np 128 /bin/true  > reportBindings.txt 2>&1

