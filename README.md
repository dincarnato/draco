![DRACO logo](http://www.incarnatolab.com/images/software/draco.png)
<br />
## Introduction

RNA structures are highly diverse and dynamic. Inside the cell, the same RNA is present in multiple identical copies, and not all of them will fold into the same structure. Rather, many alternative structures (or conformations) for the same RNA coexist in a dynamic equilibrium. In addition to that, each conformation is not static over time, but it can interconvert between any other potential conformation for that RNA. This structurally heterogeneous set of RNAs is commonly referred to as an *ensemble*.

<u>D</u>econvolution of <u>R</u>NA <u>A</u>lternative <u>CO</u>nformations (or __DRACO__), is a method that, by using a combination of spectral deconvolution and fuzzy clustering, can deconvolute the number of coexisting RNA structural conformations in a biological sample, starting from mutational profiling (MaP) data, as well as reconstruct their relative stoichiometries and individual reactivity profiles.

In order to use DRACO, BAM files from mutational profiling experiments must be pre-processed into  mutation map (MM) files using the [__RNA Framework__](https://github.com/dincarnato/RNAframework). 

For support requests, please post your questions to: <https://github.com/dincarnato/draco/issues>

For a complete documentation, please refer to: <https://draco-docs.readthedocs.io/en/latest/>


## Author(s)

Edoardo Morandi (emorandi[at]rnaframework.com) [Main developer]<br/>
Danny Incarnato (dincarnato[at]rnaframework.com)<br/>


## Reference

Morandi *et al*., (2021) Genome-scale deconvolution of RNA structure ensembles (PMID: [33619392](https://pubmed.ncbi.nlm.nih.gov/33619392/)).


## License

This program is free software, and can be redistribute and/or modified under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

Please see <http://www.gnu.org/licenses/> for more information.


## Prerequisites

- Linux system
- GCC v14 or higher (<https://gcc.gnu.org/gcc-14/>), or
  <br/>Clang v10 or greater (<https://releases.llvm.org/download.html>)
- CMake v3.8 or higher (<https://cmake.org/download/>)
- OpenBLAS (<https://www.openblas.net>)
- Armadillo v9.850.1 or greater (<http://arma.sourceforge.net>)
- Boost v1.66 or greater (<https://www.boost.org>)
- dlib v19.4.0 or greater (<http://dlib.net>)
- Intel oneAPI Thread Building Blocks (<https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onetbb.html>)

To compile the `simulate_mm` utility, Rust and Cargo are also needed (<https://doc.rust-lang.org/cargo/getting-started/installation.html>).


## Installation

Clone the DRACO git repository:

```bash
git clone https://github.com/dincarnato/draco
```
This will create a "draco" folder.<br/>
To compile DRACO:

```bash
cd draco
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DLINK_TIME_OPTIMIZATIONS=ON -DNATIVE_BUILD=ON -DARMA_NO_WRAPPER=ON
make -jN
```
where *N* is the number of processors on your computer. The parameter `CMAKE_INSTALL_PREFIX` can be used to set the installation directory. For example, to install DRACO in `/usr/local/bin`, just set it to `/usr/local`.<br/>
The `draco` executable will be located under `build/src/`.
<br/>
To compile the `simulate_mm` utility:

```bash
cd extra/simulate_mm
cargo build --release
```
The `simulate_mm` executable will be located under `target/release/`.


## Avoiding BLAS-related issues

In some cases, we observed issues related to the internal multithreading used by OpenBLAS. To avoid these and to make sure that only DRACO's multithreading is used, set the following environment variables:

```bash
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export OMP_NUM_THREADS=1
```


## Testing your installation

Under `examples/` it is possible to find three sample MM files, each one containing simulated data for a single transcript, forming one, two, or three coexisting conformations.<br/>
To analyze these samples, simply run:

```bash
draco --mm 1conf.mm        # 1 conformation
draco --mm 2confs.mm       # 2 conformations
draco --mm 3confs.mm       # 3 conformations
```
