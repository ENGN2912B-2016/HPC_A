#!/bin/bash

# Request an hour of runtime:
#SBATCH --time=20:00:00

# Default resources are 1 core with 2.8GB of memory.

#SBATCH --mem=48G

# Specify a job name:
#SBATCH -J ternary_serial_equal_ie

# Specify an output file
#SBATCH -o ternary_serial_equal_ie.out
#SBATCH -e ternary_serial_equal_ie.out

# Run a command

mkdir output
./ternary_serial 4.0 4.0 4.0

