# GROMACS on Android via Termux: MD Simulation Setup

- This project contains essential files for running Molecular Dynamics (MD) simulations using **GROMACS** on **Android** via **Termux**.

- This repository provides a minimal, tested setup for compiling **GROMACS 2024** on Termux.

## Prerequisites

This setup was **tested on Termux v0.119.0-beta.3** (downloaded from [GitHub Releases](https://github.com/termux/termux-app/releases)).  
Ensure you're using a compatible version to avoid build issues.

Install dependencies:

```bash
pkg update && pkg upgrade
pkg install root-repo x11-repo
pkg install clang cmake build-essential wget which binutils git
````

## GROMACS 2024 Compilation on Termux

Follow the steps below inside your Termux shell:

Download GROMACS source
```bash
wget https://ftp.gromacs.org/gromacs/gromacs-2024.tar.gz
tar xfz gromacs-2024.tar.gz
cd gromacs-2024
```

Fix compatibility issues
```bash
sed -i '/pthread_cancel/d' src/external/thread_mpi/src/pthreads.cpp
sed -i '/getlogin_r/d' src/gromacs/utility/sysinfo.cpp
``` 

Build GROMACS
```bash
mkdir build && cd build
cmake .. \
  -DREGRESSIONTEST_DOWNLOAD=ON \
  -DGMX_THREAD_MPI=OFF \
  -DGMX_SIMD=ARM_NEON_ASIMD \
  -DGMX_BUILD_OWN_FFTW=ON \
  -DCMAKE_INSTALL_PREFIX=$HOME/gromacs

make -j$(nproc)
make check
make install
```

Add GROMACS to your PATH
```bash
echo 'export PATH=/data/data/com.termux/files/home/gromacs/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

## Example Usage

```bash
git clone https://github.com/kssrikar4/Molecular-Dynamics-with-GROMACS.git
cd Molecular-Dynamics-with-GROMACS
#Change to your target protein
wget https://files.rcsb.org/download/3gv3.pdb -O protein.pdb

gmx pdb2gmx -f protein.pdb -o processed.gro -water tip3p -ignh
gmx editconf -f processed.gro -o newbox.gro -c -d 1.0 -bt cubic
gmx solvate -cp newbox.gro -cs spc216.gro -o solvated.gro -p topol.top
gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o solvated_ions.gro -p topol.top -pname NA -nname CL -neutral

gmx grompp -f em.mdp -c solvated_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt
gmx grompp -f md.mdp -c npt.gro -p topol.top -o md.tpr
gmx mdrun -deffnm md
```

## Operational Screenshots 

### GROMACS initial execution displaying synopsis and options within Termux on Android
![termux1](https://github.com/kssrikar4/Molecular-Dynamics-with-GROMACS/blob/350a5b378f9888eb5540b9697e0f3028cdd4f5ce/plots/Termux-1.jpg) 

### GROMACS execution in Termux, showing simulation parameters and real-time progress for a molecular dynamics run.
![termux2](https://github.com/kssrikar4/Molecular-Dynamics-with-GROMACS/blob/350a5b378f9888eb5540b9697e0f3028cdd4f5ce/plots/Termux-2.jpg)

### Monitoring GROMACS process using `top` in Termux, displaying CPU and memory utilization.
![termuxtop](https://github.com/kssrikar4/Molecular-Dynamics-with-GROMACS/blob/350a5b378f9888eb5540b9697e0f3028cdd4f5ce/plots/Termux-top.jpg)


## Disclaimer

* **Performance**: Mobile CPUs are limited so use short time steps or fewer atoms for testing.
* **Battery & Heat**: Monitor device temperature; use power banks for sustained runs.
* **Offload**: Use this as a setup for quick preprocessing, then migrate `.tpr` files to HPC clusters.

For those who want **high-performance science on pocket-sized hardware**.
