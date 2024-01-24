PSI -- Pan-genome Seed Index
============================

This is an implementation of the method introduced in:

> Ghaffaari, A. & Marschall, T.
> Fully-sensitive seed finding in sequence graphs using a hybrid index 
> _Bioinformatics_, **2019**, _35_, i81-i89

How to Install
--------------

### Dependencies

#### Minimum requirements

PSI uses many features introduced in C++17. In order to work with it, you need a
C++17 compliant compiler first. It has been tested for these compilers:

| Compiler   | Version                  |
| ---------- | ------------------------ |
| GCC        | 8.3, 9.2, 10.2.0, 12.2.0 |
| AppleClang | 14.0.0, 16.0.1           |

#### Required dependencies

PSI library is a header-only library. However, it depends on a few external
libraries to work. Both `psikt` seed finder and PSI library need following
libraries to be installed before building:

- Zlib
- BZip2
- Protobuf (tested with <= 3.21.12)

On Ubuntu they can be installed by:

``` bash
$ sudo apt install zlib1g-dev libbz2-dev libprotobuf-dev protobuf-compiler
```

---

**NOTE**
There is a bug in Protobuf 3.20.x causing
[a linking issue](https://github.com/protocolbuffers/protobuf/issues/9947) while
compiling in debug mode (i.e. when `CMAKE_BUILD_TYPE` matches `Debug`).
Otherwise, in Release mode, it should be fine.

---

Configuring the source code and/or building auxiliary tools requires:

- CMake >=3.10

#### Bundled dependencies

There are a few more dependencies that is bundled with PSI; i.e. PSI can install them
for you. If you want PSI to to do so, specify `-DUSE_BUNDLED_ALL=on` when running cmake
(explained later). These dependencies are:

- SeqAn, version 2.4.0
- GUM, version 1.2.2
- kseq++, version 1.0.0
- PairG, version 0.2.1

PSI build script will skip those dependencies that are already installed on your system.
Besides, installing each dependency can be individually enabled or disabled by setting
`-DUSE_BUNDLED_<library>` to `on` or `off` respectively.

---

**NOTE**
If the bundled version of SeqAn did not work for any reason, install it manually. Make
sure that `seqan-2.pc` configuration file has been installed too. On Debian-based
operating systems, the easiest way to install SeqAn might be using the package manager:

```bash
$ sudo apt install libseqan2-dev
```

---

### Installation

To install PSI (recommended):

```bash
$ git clone https://github.com/cartoonist/psi.git
$ mkdir psi/build
$ cd psi/build
$ cmake -DCMAKE_BUILD_TYPE=Release -DKokkos_ENABLE_OPENMP=TRUE -DUSE_BUNDLED_ALL=on ..
$ make
$ sudo make install
```

Use PSI
-------

There are two ways to include PSI into your project:

### As a git submodule

If you are using CMake, the library can be added to your project as a git
submodule. Then, calling `add_subdirectory(path/to/psi)` in the corresponding
`CMakeLists.txt` file exports `psi::psi` target. The exported target defines
include directories, transitive dependencies, and compiler flags required for
building the project.

### TODO: As an external dependency

Install the PSI as an external dependency. If you are using CMake `find_package`
will import `psi::psi` target to which can be used to link your project.

If you are not using CMake, the same information can also be retrieved using
`pkg-config`. Or just take a look at `<INSTALL_PREFIX>/lib/pkgconfig/psi.pc`.

Development
-----------

In order to build the tests or auxiliary tools, just turn on `BUILD_TESTING` or
`BUILD_PSI_AUX_TOOLS` options when running `cmake`:

```bash
$ cmake -DCMAKE_BUILD_TYPE=Debug -DKokkos_ENABLE_OPENMP=TRUE -DBUILD_TESTING=on -DBUILD_PSI_AUX_TOOLS=on ..
```
