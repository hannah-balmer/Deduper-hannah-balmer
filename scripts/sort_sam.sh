#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp

set -e # Stop task when error is encountered

/usr/bin/time -v pixi run samtools sort /projects/bgmp/shared/deduper/C1_SE_uniqAlign.sam -o /projects/bgmp/hbalmer/bioinfo/Bi624/Deduper-hannah-balmer/sorted_insam.sam