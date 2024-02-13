#!/bin/bash

#should probably make this better...
cd graph/ || {
    exit 1
}
mkdir build || {
    exit 1
}
cmake -S . -B ./build || {
    exit 1
}
cd build/ || {
    exit 1
}
make -j 4 || {
    exit 1
}
