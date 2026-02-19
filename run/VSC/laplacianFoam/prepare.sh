#!/bin/bash
#SBATCH -A p......
#SBATCH -J OF_prep
#SBATCH -t 00:10:00
#SBATCH --partition=zen3_0512
#SBATCH --qos=zen3_0512
#SBATCH --ntasks=1
#SBATCH --mem=60G

source $HOME/OpenFOAM/OpenFOAM-12/etc/bashrc

blockMesh > log.blockMesh

decomposePar -cellProc -force > log.decomposePar
