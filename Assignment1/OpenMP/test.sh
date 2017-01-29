#!/bin/bash

for size in 1024 2048 4096 8192; do
    for threads in 1 2 4 8 16; do
        echo "Threads = $threads size = $size"
        make clean
        make test NUM_THREADS=$threads MATRIX_SIZE=$size
    done
done
