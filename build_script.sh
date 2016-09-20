#!/bin/bash

## Note: the environment variable MITGCMROOT needs to be set outside this script
cd build
$MITGCMROOT/tools/genmake2 -mpi -rootdir=$MITGCMROOT -mods=../code_new -of=../code_new/linux_ia64_ifort+mpi_yeti
make depend
make -j4
