# Topology-Aware Communication Optimization for OpenFOAM

This repository provides an **OpenFOAM-12 utility** for topology-aware communication optimization in parallel CFD simulations on NUMA-based HPC systems.

The tool improves communication performance by:

- Optimizing **MPI rank-to-core placement** based on the communication graph induced by the mesh decomposition
- Introducing a **hierarchical global reduction schedule** aligned with the hardware topology

This implementation accompanies the paper:

> *Topology-Aware Communication Optimization for CFD Simulations in OpenFOAM*

---

## Motivation

In OpenFOAM, the finite-volume discretization leads to large sparse systems of linear equations that are solved using iterative methods (e.g., CG, GAMG). These solvers dominate runtime and require substantial:

- **Local communication** (halo exchange during sparse matrix–vector products)
- **Global communication** (collective reductions such as scalar products and norms)

The default MPI rank assignment and reduction schedules in OpenFOAM are largely topology-agnostic. On NUMA-based cluster architectures, this can lead to communication across slow inter-socket or inter-node links, even when faster local alternatives exist.

This utility analyzes the processor patch structure of a given decomposition and generates:

- An optimized MPI **rank file**
- A topology-aware **global reduction schedule**

---

## Requirements

- OpenFOAM-12
- Compatible MPI implementation (e.g., OpenMPI 4.x)
- Linux-based HPC system (SLURM recommended)
- SCOTCH (from OpenFOAM `ThirdParty-12`)
- Homogeneous compute nodes (assumed)

---

## Installation

### 1. Copy the Utility

Create the directory 

```
$WM_PROJECT_USER_DIR/applications
```
Copy the `optimizeCommPar` directory into it.



### 2. Compile the Utility

```bash
./Allwmake
```

### 3. Enable Custom Global Communication

Replace the following OpenFOAM source files:

```
src/OpenFOAM/db/IOstreams/Pstreams/UPstream.H
src/OpenFOAM/db/IOstreams/Pstreams/UPstream.C
src/OpenFOAM/db/IOstreams/Pstreams/PstreamReduceOps.H
```

Then recompile OpenFOAM:

```bash
wmake
```

⚠️ Note: Modifying OpenFOAM core files is required to enable custom reduction schedules.

---

## Case Setup

### Local Communication Optimization

Add the following sub-dictionary to:

```
system/decomposeParDict
```

Example:
```cpp
optimizeCommCoeffs
{
    procPerNode     24;
    numaDomains     2;
    numaNodes       4;

    version         0;

    domainList      (0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1);
    numaList        (0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3 0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3);
}
```

### Parameter Description

- `procPerNode` – Number of cores per compute node  
- `numaDomains` – Number of NUMA domains per compute node  
- `numaNodes` – Number of NUMA nodes per NUMA domain  

The number of cores per NUMA node is therefore

```
procPerNode / (numaDomains * numaNodes)
```

The arrays `domainList` and `numaList` must accurately reflect the actual hardware topology. NUMA domains are indexed from `0` to `numaDomains - 1`. Within each NUMA domain, NUMA nodes are indexed from `0` to `numaNodes - 1`.

These lists must represent the exact mapping of cores to NUMA domains and NUMA nodes. The required hardware information can be obtained, for example, using:

```
lscpu
numactl --hardware
```

### Optimization Modes

| Version | Description |
|----------|------------|
| `0` | Local + Global communication optimization |
| `1` | Random rank file (testing only) |
| `2` | Global communication optimization only |

---

## Global Communication Optimization

Add to:

```
system/controlDict
```

```cpp
OptimisationSwitches
{
    commSchedule 3;
}
```

---

## Running Cases

Prepare the case as usual:

```bash
blockMesh
decomposePar
```

---

### Local + Global Optimization

For local and global communication optimization, first execute the utility in the same SLURM job script. 
Then start the OpenFOAM solver with the generated rank file using `-rf ./constant/rankFile`, for example:

```bash
mpirun -np 24 optimizeCommPar -parallel > log.optimizeCommPar
mpirun -rf ./constant/rankFile -np 24 foamRun -parallel > log.foamRun
```

---

### Global Optimization Only

1. Determine the MPI rank-to-core mapping using `reportBindings.sh`
2. Place `reportBindings.txt` in the case directory
3. Run:

```bash
mpirun -np 24 optimizeCommPar -parallel > log.optimizeCommPar
mpirun -bind-to core -np 24 foamRun -parallel > log.foamRun
```

The SLURM job headers must be adapted to your HPC system.

---

## Test Cases

The `run/` directory contains the four test cases presented in the paper, configured for the respective HPC systems.
It also includes the corresponding `prepare.sh` setup script and SLURM job scripts for the different optimization modes.

The provided case configurations are set up for execution on two compute nodes. 
For different node counts, only the `decomposeParDict` and the corresponding SLURM job scripts need to be adapted accordingly.

---

## Important Notes

- Correct hardware topology specification is essential.
- The implementation assumes **homogeneous compute nodes**.
- Incorrect NUMA mappings may degrade performance.
- The optimization is static and assumes a fixed domain decomposition.

---

## Disclaimer

This software modifies internal OpenFOAM communication components.  
Use at your own risk and validate results carefully.
