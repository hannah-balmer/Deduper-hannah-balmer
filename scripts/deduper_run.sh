#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp

set -e # Stop task when error is encountered


/usr/bin/time -v pixi run ./balmer_deduper.py \
    -u ../files/STL96.txt \
    -f ../files/sorted_insam.sam \
    -o ../files/sorted_outsam.sam

"""
/usr/bin/time -v pixi run ./balmer_deduper.py \
    -u ../files/STL96.txt \
    -f ../unit_tests/input.sam \
    -o ../unit_tests/test_output.sam
"""

