# zfp_example

This repository shows some basic usage of zfp:

https://computing.llnl.gov/projects/floating-point-compression

https://github.com/LLNL/zfp

Compilation
-----------
Before compiling, set `ZFP_DIR` and `ZFP_LIB` environment variables.
On most UNIX-based systems a default compilation can be performed using the
following commands from the top-level of the source tree:

    mkdir build
    cd build
    cmake ..
    make
