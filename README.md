# DRACO: Deconvolution of RNA Alternative COnformations

## Requirements

* A linux-based distribution
* GCC 9 or higher (GCC 10 is better because GCC 9 gives a lot of false positive warnings), or Clang 9 or higher
* [Armadillo](http://arma.sourceforge.net/) library
* [Dlib](http://dlib.net/) library
* [Intel TBB](https://software.intel.com/content/www/us/en/develop/tools/threading-building-blocks.html)

It is highly recommended to use the package manager of your distribution to install the required dependencies.

## How to compile

The project uses [CMake](https://cmake.org/) as building tool.

After cloning the repository you can do the following:

```bash
cd draco
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DLINK_TIME_OPTIMIZATIONS=ON -DNATIVE_BUILD=ON -DARMA_NO_WRAPPER=ON
make -jN
```

Where `N` is generally not above the number of the cores on your machine.

It is also possible to use `ccmake` in order to customize some building variables.

The `draco` executable can be found inside `build/src/`, and it is also possible to install it in `/usr/local/bin` (or where the `${CMAKE_INSTALL_PREFIX}/bin` points to) using `make install`.

## Running DRACO

*DRACO* requires a DMS-MaPseq data experiment that must be processed using the [RNA Framework](https://www.rnaframework.com/). A BAM files of mapped reads must be given to the `rf-count` tool to produce MM files, which can be analyzed using *DRACO*.
