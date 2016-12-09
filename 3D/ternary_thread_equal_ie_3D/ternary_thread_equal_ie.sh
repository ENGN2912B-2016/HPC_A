#!/bin/bash

# Request an hour of runtime:
#SBATCH --time=48:00:00

# Default resources are 1 core with 2.8GB of memory per core.

# Use more cores (8):
#SBATCH -c 12
#SBATCH --mem=48g

# Specify a job name:
#SBATCH -J ternary_thread_equal_ie_3D

# Specify an output file
#SBATCH -o ternary_thread_equal_ie_3D.out
#SBATCH -e ternary_thread_equal_ie_3D.out

# Run a command

mkdir output
./ternary_thread_3D 4.0 4.0 4.0 -threads 12

