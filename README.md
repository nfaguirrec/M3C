# M3C

M3C is a general-purpose code, though its primary targets are those fragmentation processes that take place through non- radiative transitions, where their fragments remain together enough time such that the excess of energy becomes randomly distributed over all internal degrees of freedom. M3C allows studying fragmentation processes in a large variety of systems irrespective of their composition and the nature of the chemical bonding between their elementary constituents.

In essence, the M3C consists in partitioning the mass, charge, energy, and momentum (linear and angular) of an excited molecular system (which are conserved in the microcanonical approach) among all accessible fragmentation channels with probabilities governed by considerations of maximum entropy. The key aspect of this methodology is that it provides a random way to move in phase space until a region of maximum entropy is reached, where the physical observables are computed. Here, only energy and entropy are responsible for the fate of the system.

The simplest calculation requires quantities such as geometries, electronic energies, and harmonic frequencies of all molecules included in the fragmentation model. This approach has already successfully been applied to the fragmentation of neutral and singly-charged carbon clusters, hydrogenated and nitrogenated carbon clusters, and small and medium-sized molecules in different physical scenarios or ionizing radiation.

## Requirements:

Required:

- **[Intel® Fortran Compiler](https://software.intel.com/en-us/fortran-compilers)** (version >= 14.0.3). M3C has not been tested with any other compiler.

- **[SciFT (Scientific Fortran Tools)](https://github.com/nfaguirrec/scift)**<br>
  The Scientific Fortran Tools (SciFT) is a numerical library for fortran programmers.

Optional:

- **[libmsym](https://github.com/mcodev31/libmsym)**<br>
  libmsym is a C library dealing with point group symmetry in molecules.

## Compiling M3C

Download the .zip file from this page and extract the files,
```
$ unzip M3C-master.zip 
Archive:  M3C-master.zip
9f1572142803f97705c5db16b2d018a9c853a658
   creating: M3C-master/
  inflating: M3C-master/LICENSE      
  inflating: M3C-master/LICENSE.jmol
...

$ mv M3C-master M3C
```
or clone the repository using git
```
$ git clone https://github.com/nfaguirrec/M3C.git
```
The following should be the content of the M3C-store directory if previous steps were successful:
```
$ cd M3C
$ ls
doc   doxyfile  LICENSE.jmol     LICENSE.scift  Makefile   src        utils
docs  LICENSE   LICENSE.libmsym  M3Cvars.sh     README.md  templates  VERSION
```

Enter in the M3C directory (cd M3C) and modify the Makefile file (src/Makefile). In particular, choose the right path to the scift library (-I${HOME}/Software/scift/src and -L${HOME}/Software/scift/src).

To build the code just type make inside the main directory as follows:
```
$ make
cd src; make; cd ..
make[1]: Entering directory '/scratch/nestor/M3C/src'
Building dependencies for AzizSlamanPotential.f90 ... OK
Building dependencies for Fragment.f90 ... OK
Building dependencies for FragmentsDB.f90 ... OK
...
Building MarkovChain.f90 (0:03.40)
Building NNLS.f90 (0:00.08)
Building M3CBR.f90 (0:01.53)
```

## Installing M3C

The basic environmental variables that M3C needs can be loaded just adding the following command anywhere in the ~/.bashrc file:

```
source <PATH_TO_M3C>/M3Cvars.sh
```

M3C is also able to obtain data from electronic structure calculations by interfacing with some standard quantum chemistry programs. To enable this option the following variables must be specified:

```
# GAMESS configuration
export M3C_GAMESS_HOME=<PATH_TO_GAMESS_INSTALLATION>
export M3C_GAMESS_SCRATCH=/scratch/$USER/gamess

# GAUSSIAN configuration
export M3C_GAUSSIAN_HOME=<PATH_TO_GAUSSIAN_INSTALLATION>
export M3C_GAUSSIAN_SCRATCH=/scratch/$USER/gaussian
```

# Authors

* Nestor F. Aguirre ( nfaguirrec@gmail.com )
* Sergio Díaz-Tendero ( sergio.diaztendero@uam.es )
* Paul-Antoine Hervieux ( paul-antoine.hervieux@ipcms.unistra.fr )
* Manuel Alcamí ( manuel.alcami@uam.es )
* Fernando Martín ( fernando.martin@uam.es )

## Citing

1. Fully versus Constrained Statistical Fragmentation of Carbon Clusters and their Heteronuclear Derivatives.<br/>
   N. F. Aguirre, S. Díaz-Tendero, T. IdBarkach, M. Chabot, K. Béroff, M. Alcamí, and F. Martín.<br/>
   [*J. Chem. Phys.* **150**, 144301 (2019)](https://aip.scitation.org/doi/10.1063/1.5083864)

2. Furan Fragmentation in the Gas Phase: New Insights from Statistical and Molecular Dynamics Calculations.<br/>
   E. Erdmann, M. Łabuda, N. F. Aguirre, S. Díaz-Tendero, and M. Alcamí.<br/>
   [*J. Phys. Chem. A* **122**, 4153-4166 (2018)](http://pubs.acs.org/doi/10.1021/acs.jpca.8b00881)

3. Semiempirical breakdown curves of C2N(+) and C3N(+) molecules; application to products branching ratios predictions of physical and chemical processes involving these adducts.<br/>
   T. IdBarkach, T. Mahajan, M. Chabot, K. Béroff, N.F. Aguirre, S. Diaz-Tendero, T. Launoy, A. Le Padelle, L. Perrot, M.A. Bonnin, K.C. Le, F. Geslin, N. de Séréville, F. Hammache, A. Jallat, A. Meyer, E. Charon, T. Pino, T. Hamelin, and V. Wakelam<br/>
   [*Mol. Astrophysics* **12**, 25-32 (2018)](https://www.sciencedirect.com/science/article/pii/S2405675818300125)

5. M3C: A Computational Approach To Describe Statistical Fragmentation of Excited Molecules and Clusters.<br/>
   N. F. Aguirre, S. Díaz-Tendero, P.-A. Hervieux, M. Alcamí, and F. Martín.<br/>
   [*J. Chem. Theory Comput.* **13**, 992-1009 (2017)](http://pubs.acs.org/doi/pdf/10.1021/acs.jctc.6b00984)

5. Structure, Ionization, and Fragmentation of Neutral and Positively Charged Hydrogenated Carbon Clusters:
   C<sub>n</sub>H<sub>m</sub><sup>q+</sup> (n=1-5, m=1-4, q=0-3).<br/>
   J. P. Sánchez, N. F. Aguirre, S. Díaz-Tendero, F. Martín, and M. Alcamí.<br/>
   [*J. Phys. Chem. A* **120**, 588-605 (2016)](http://pubs.acs.org/doi/abs/10.1021/acs.jpca.5b10143)
