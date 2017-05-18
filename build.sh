#!/usr/bin/env sh

mkdir build
cd build
cmake ..
make -j
make test
cd ..
