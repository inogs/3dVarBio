# Parallel implementation of a data assimilation scheme for operational oceanography: the case of the MedBFM model system

By Teruzzi, A., Di Cerbo P., Cossarini, G., Pascolo, E., Salon, S.

### Abstract

The MedBFM model system provides forecasts and reanalysis of the Mediter- ranean Sea biogeochemistry for the European Copernicus Services. The system integrates model and observations through a 3D variational assimilation scheme, whose performance capitalizes on present HPC systems and ensures compli- ance to the service requirements. Domain decomposition with message passing paradigm was implemented to parallelize the assimilation code and maximize performance and scalability. In particular, an innovative approach to the opti- mal parallelization of the spatial filtering algorithm of the assimilation scheme with a dynamical sliced decomposition was developed. Moreover, the efficient parallel solver of the PETSc/TAO library was adopted for optimizing the cost function minimization on which the variational assimilation is based. Consider- ing the complex domain decomposition of the Mediterranean Sea and the high variability of the data input, a modern and scalable software application for variational assimilation schemes that exploits the performance of modern HPC architectures on whole node was developed.

## 3DVarBio instruction guide

### Download

The source files of the 3DVarBio code can be downloaded from the following git repository: <https://github.com/inogs/3dVarBio.git>.

The version to be used is tagged “v1.0” and can be downloaded with the following git instructions:

```
git clone git@github.com:inogs/3dVarBio.git 
```

```
cd 3dVarBio
```

```
git checkout tags/v1.0
```

### Code structure

The `main.f90` subroutine of the 3DVarBio code calls the subroutine `biovar.f90`, which then calls all the subroutines for the execution of the code.
The `namelists/` directory contains an example of a namelist file suitable for this code version and to configure the launch of a test run. The main settings of the namelist file and the values used in the test run provided are the following:

* `ctl_tol`
   Absolute stopping criteria (not used in TaoMinimizer)

* `ctl_per`
   Relative stopping criteria (the minimization ends when the gradient of the cost function is less then ctl_per of the gradient at the first step)

* `Neof` Number of EOF profiles (Vv operator)

* `Nreg` Number of regions with different EOF preofiles (Vv operator)

* `rcf_ntr` Number of iterations of the recursive filter (Vh operator)

* `rcf_L` Horizontal correlation radius in meters (Vh operator)

* `rcf_efc` Extension factor for coasts (Vh operator)

* `bio_assim` Updates of the assimilation for all the involved biogeochemical variables (Vb operator): =1 calculated internally to 3DVarBio with output made by restart files for the biogeochemical model; =0 not calculated  and 3DVarBio provides updates only for total chlorophyll without applying conditions to check updates consistency (output made by a single corr.nc file)

* `Nphyto` Number of phytoplankton functional types

* `chldep` Minimum water depth (in meters) of areas of assimilation of satellite observations

* `ncmp` Number of components of phytoplankton functional types

* `ApplyConditions` In the evaluation of the updates for the biogeochemical variables 3DVarBio (Vb operator) applies conditions (ApplyCponditions=.true.) to check the consistency of the updates (used only if bio_assim=1)

* `sat_obs` Assimilation of satellite observations (sat_obs=1)

* `argo` Assimilation of argo float observations (argo=1)

* `uniformL` Use of non-uniform correlation radius (uniformL=1; Vh operator)

* `anisL` Use of anisotropic correlation radius (anisL=1; Vh operator)

* `verbose` Set verbose output (verbose=1)

### Compilation environment

The compilation was tested with an intel compiler and using the following modules on the HPC cluster PICO  (66 nodes x Intel Xeon E5 2670 v2 @2.5Ghz 128 GB, Mellanox Infinitband FDR) located at CINECA supercomputing center (Bologna, Italy; <http://www.hpc.cineca.it/>):



```
profile/advanced
```

```
autoload/0.1
```

```
intel/pe-xe-2016—binary
```

```
intelmpi/5.1.3—binary
```

```
petsc/3.7.2--intelmpi--5.1.3—binary
```

```
netcdf/4.4.0-parallel--intelmpi--5.1.3—binary
```

```
netcdff/4.4.4-parallel--intelmpi--5.1.3—binary
```

```
pnetcdf/1.7.0--intelmpi--5.1.3—binary
```

For the compilation it is necessary to prepare the file compiler.inc. This file contains all the information about the compilers, the libraries and the flags adopted during the compilation. In our configuration we used x86_64.LINUX.intel.inc that is available in the repository. Other templates are available within the repository (x86_64.LINUX.*). We suggest to copy one of the templates (e.g., x86_64.LINUX.intel.inc) into compiler.inc, and then to modify it accordingly to the user configuration.
A copy of an execution test can be found in the following repository <https://github.com/inogs/3DVarbio_test.git>.

### Test run: details on configuration and launch

The test was carried out on single node using 20 cores. The files contained in the test directory are provided herafter.

#### Input files

* `eofs.nc` Containing a set of EOFs which accounts for vertical covariance (Vv)

   `eva` EOF eigenvalues

   `evc` EOF eigenvector

* `grid1.nc` Grid information 
   
   `dep` grid levels

   `dx` horizontal resolution along longitude
   
   `dy` horizontal resolution along latitude
   
   `lat` latitude
   
   `lon` longitude
   
   `regs` identification index of regions with homogenous EOFs
   
   `tmsk` mask of the grid

* `chl_mis.nc` Misfit between satellite and model chlorophyll and observation error
   
   `misfchl` misfit y-Hx
   
   `errchl` observation error R

* `gradsal.nc` Provide the parameters for the use of anisotropic correlation radius (Vh)
   
   `kx_n` factor for the radius along longitude

   `ky_n` factor along latitude

* `chl_rad_corr.nc` Non homogenous correlation radius

   `radius`

* `DA__FREQ_1/RST.20130101-120000.P??.nc` 17 files for each component of the 4 phytoplankton functional types for all the vertical levels of the biogeochemical model. The name of the files refers to the name of the contained variable Pjk. Where, j is the number of the phytoplankton functional type (ranging from 1 to 4);
k identifies the internal component (c: carbon; l: chlorophyll; n: nitrogen; p: phosphorous; s: silicon only for phytoplankton P1, i.e., diatoms)

* `var_3d_nml`

#### Output files

* `RESTARTS/RST.20130101-120000.P??.nc` Updates provided by 3DvarBio for the 17 biogeochemical variables that describes the phytoplankton functional types. The name of the file is made using Pjk, while the variable name in each file is TRNPjk:
j is the number of the phytoplankton functional type (ranging from 1 to 4);
k identifies the internal component (c: carbon; l: chlorophyll; n: nitrogen; p: phosphorous; s: silicon only for phytoplankton P1, i.e., diatoms)


* `RESTARTS/RST.20130101-12:00:00.P??.nc` Symbolic link to the previous files needed to correctly couple 3DVarBio with the biogeochemical model (OGST-BFM)

* `BioVar.diagnostics` File containing diagnostics of the 3DVarBio run

#### Program file
* `var_3d`

The code has been executed using the following command:

```
mpirun –np 20 ./var_3d
```




