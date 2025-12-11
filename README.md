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


## Installation

Clone the DRACO git repository:

```bash
git clone https://github.com/dincarnato/draco
```
This will create a "draco" folder.

There are two ways you can compile and use DRACO.

### Using Docker (recommended)

You can use [docker](https://www.docker.com/) to build an image with all the necessary pre-requisites to compile and run DRACO.

#### Building the image

```shell
cd draco
docker build --tag draco .
```

#### Running DRACO

To run DRACO using the docker image, it is necessary to expose the input and output directories to the docker container that is being spawned.

For instance, if your `.mm` files are located in the `/my/mm/files` folder, docker needs to be invoked as it follows:

```shell
docker run --init --rm -u $(id -u ${USER}):$(id -g ${USER}) -v /my/mm/files:/data draco --mm rna.mm <other draco parameters>...
```

where:

- `--init` forwards the signals to the process (i.e., stopping using CTRL+C)
- `--rm` cleans up the container after its use
- `-u $(id -u ${USER}):$(id -g ${USER})` runs DRACO with the specified user and group id
- `-v /my/mm/files:/data` is to mount the directory `/my/mm/files` to `/data` inside the container. DRACO runs inside this directory, therefore all the paths should be relative and within the directory `/my/mm/files` (see the next parameter)
- `--mm rna.mm` is the MM file to be analyzed, which must be located inside the `/my/mm/files` directory

DRACO's output files will be placed inside the `/my/mm/files` directory.

### Compiling and running DRACO without Docker

#### Prerequisites

- Linux system
- GCC v15 or higher (<https://gcc.gnu.org/gcc-15/>), or
  <br/>Clang v19 or greater (<https://releases.llvm.org/download.html>)
- CMake v3.8 or higher (<https://cmake.org/download/>)
- OpenBLAS (<https://www.openblas.net>)
- Armadillo v9.850.1 or greater (<http://arma.sourceforge.net>)
- Boost v1.66 or greater (<https://www.boost.org>)
- dlib v19.4.0 or greater (<http://dlib.net>)
- Intel oneAPI Thread Building Blocks (<https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onetbb.html>)

#### Compilation

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


### Compiling and running the `simulate_mm` utility

#### Prerequisites

To compile the `simulate_mm` utility, Rust and Cargo are needed (<https://doc.rust-lang.org/cargo/getting-started/installation.html>).

#### Compilation

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
To analyze these samples, if DRACO is compiled within Docker simply run:

```bash
docker run --init --rm -u $(id -u ${USER}):$(id -g ${USER}) -v examples:/data draco --mm 1conf.mm        # 1 conformation
docker run --init --rm -u $(id -u ${USER}):$(id -g ${USER}) -v examples:/data draco --mm 2confs.mm       # 2 conformations
docker run --init --rm -u $(id -u ${USER}):$(id -g ${USER}) -v examples:/data draco --mm 3confs.mm       # 3 conformations
```

If DRACO is compiled as a stand-alone executable, run:

```bash
./build/src/draco --mm examples/1conf.mm        # 1 conformation
./build/src/draco --mm examples/2confs.mm       # 2 conformations
./build/src/draco --mm examples/3confs.mm       # 3 conformations
```

## Known issues

### Clang &ge; 19 and dlib

On Ubuntu 25.04, if Clang is used to compile DRACO and LibC++ is used (with `-DUSE_LIBCXX=ON`), a compilation error is triggered, which appears to be caused by dlib.

The error message looks like the following:
```
/usr/include/dlib/matrix/../serialize.h:1963:14: note: in instantiation of member function 'std::basic_string<unsigned int>::resize' requested here
 1963 |         item.resize(size);
      |              ^
```

This has been reported in [a closed dlib issue](https://github.com/davisking/dlib/issues/3045), and the problem should be fixed in version 19.24.7.

Unfortunately, to date, as Ubuntu 25.04 ships with an older version of the library, we suggest to use GCC or, if possible, to manually update dlib.
