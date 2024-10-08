# ParsiMoNe - Parallel Construction of Module Networks
[![Build](https://github.com/asrivast28/ParsiMoNe/actions/workflows/main.yml/badge.svg)](https://github.com/asrivast28/ParsiMoNe/actions/workflows/main.yml)
[![Apache 2.0 License](https://img.shields.io/badge/license-Apache%20v2.0-blue.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/349758347.svg)](https://zenodo.org/badge/latestdoi/349758347)


ParsiMoNe (**Par**allel Con**s**truct**i**on of **Mo**dule **Ne**tworks) supports learning of module networks in parallel.

## Requirements
* **gcc** (with C++14 support) is used for compiling the project.  
_This project has been tested only on Linux platform, using version [10.1.0](https://gcc.gnu.org/gcc-10/changes.html)._
* **[Boost](http://boost.org/)** libraries are used for parsing the command line options, logging, and a few other purposes.  
_Tested with version [1.74.0](https://www.boost.org/users/history/version_1_74_0.html)._
* **[TRNG](https://www.numbercrunch.de/trng/)** is used for generating pseudo random numbers sequentially and in parallel.  
_Tested with version [4.22](https://github.com/rabauke/trng4/releases/tag/v4.22)._
* **[Armadillo](http://arma.sourceforge.net/)** is used for executing linear algebra operations during consensus clustering.  
_Tested with version [9.800.3](http://sourceforge.net/projects/arma/files/armadillo-9.800.3.tar.xz)._
* **[MPI](https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/mpi31-report.htm)** is used for execution in parallel.  
_Tested with [MVAPICH2 version 2.3.3](http://mvapich.cse.ohio-state.edu/static/media/mvapich/mvapich2-2.3.3-userguide.html)._
* **[CMake](http://cmake.org/)** is required for building the project.  
_Tested with version [3.29](https://cmake.org/cmake/help/v3.29/)._
* The following repositories are used as submodules:
  * **[BN Utils](https://github.com/AluruLab/bn-utils)** contains common utilities for learning in parallel and scripts for post-processing.  
  * **[mxx](https://gitlab.com/patflick/mxx)** is used as a C++ wrapper for MPI.  
  * **[C++ Utils](https://github.com/asrivast28/cpp-utils)** are used for logging and timing.  
  * **[trng4](https://github.com/rabauke/trng4/)** is used for random number generation.

## Building
After the dependencies have been installed, the project can be built as:  
<pre><code>mkdir build
cd build
cmake -DArmadillo_ROOT=${ARMA_INSTALL_LOCATION} ..
</code></pre>  
This will create an executable named `parsimone`, which can be used for constraint-based structure learning.  

#### Debug
For building the debug version of the executable, the following can be executed:
<pre><code>cmake -DCMAKE_BUILD_TYPE=Debug .. 
</code></pre>  

#### Logging
By default, logging is disabled in the release build and enabled in the debug build.
In order to change the default behavior, `LOGGING` argument can be passed to `cmake`:  
<pre><code>cmake -DENABLE_LOGGING=ON
</code></pre>
Please be aware that enabling logging will affect the performance.

#### Timing
Timing of high-level operations can be enabled by passing `-DENABLE_TIMING=ON` argument to `cmake`.

## Execution
Once the project has been built, please execute the following for more information on all the options that the executable accepts:
<pre><code>./parsimone --help
</code></pre>
For running in parallel, the following can be executed:
<pre><code> mpirun -np 8 ./parsimone ...
</code></pre>  

## Algorithms
Currently, the only supported algorithm for learning module networks is `lemontree` that corresponds to the algorithm by [Bonnet et al.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003983) originally implemented in [_Lemon-Tree_](https://github.com/erbon7/lemon-tree).

## Publication
[**Ankit Srivastava, Sriram Chockalingam, Maneesha Aluru, and Srinivas Aluru.** "Parallel Construction of Module Networks."
_In 2021 SC21: International Conference for High Performance Computing, Networking, Storage and Analysis (SC)_, IEEE Computer Society, 2021.](https://dl.acm.org/doi/10.1145/3458817.3476207)

_The experiments in the publication can be reproduced using [`EXPERIMENTS.md`](EXPERIMENTS.md)._

## Licensing
Our code is licensed under the Apache License 2.0 (see [`LICENSE`](LICENSE)).
