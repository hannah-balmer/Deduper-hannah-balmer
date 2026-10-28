# Lab_Notebook_
 
**Start date: October 24, 2026**

**Repository: /projects/bgmp/hbalmer/bioinfo/Bi624/Deduper-hannah-balmer**
 
---
 
## Project Overview
 
**Objective:**
Creation of a pipeline that takes an unsorted raw SAM file and sorts it to then remove PCR duplicates from it, outputting a new sorted and filtered SAM file and including some filtering metrics in standard output.
 
## Environment & Reproducibility
 
| Item | Details |
|---|---|
| OS | Windows 11 |
| Compute resource | Talapas |
| Package manager | Pixi |
| Environment file location | /projects/bgmp/hbalmer/bioinfo/Bi624/Deduper-hannah-balmer/pixi.toml |
| Version control repo | https://github.com/hannah-balmer/Deduper-hannah-balmer |
 


**Environment setup notes:**
```bash
pixi init
pixi add samtools
```
---
 
### Data and Overall Results

| Dataset | Source | Version | Date obtained |
|---|---|---|---|
| /projects/bgmp/shared/deduper/C1_SE_uniqAlign.sam | Talapas | Unsorted SAM file | 10/24/26 |

---
| SLURM process | Time | %CPU | Memory (kbytes)|
|---|---|---|---|
| Sorting input SAM | .25 | 96% | 934024 |
| Deduping sorted SAM | 1:09.23 | 99% | 557260 |


**Data description:**
Sequencing SAM file is single end and needs to be sorted prior to input in program. In col3 (RNAME/chromosome #), the last couple thousand reads do not have traditional chromosome names as they likely are alignemnts to the mitochondrial chromosome or other scaffolds not associated with a chromsome?

---
 
## Daily Log
 

### 2026-10-24
 
**Goal for today:**
Begin deduper script and sort input file with samtools

**Steps performed:**
1. Sorted input SAM file and do intiial data exploration with sorted file
2. Completed rough draft of deduper script 
    - Currently missing: log file output code, corrected cigar string logic, output file writing code, and testing
 
**Commands / scripts run:**
```bash
# Sorting input SAM file with sort_sam.sh:
/usr/bin/time -v pixi run samtools sort /projects/bgmp/shared/deduper/C1_SE_uniqAlign.sam -o /projects/bgmp/hbalmer/bioinfo/Bi624/Deduper-hannah-balmer/sorted_insam.sam
``` 
 
**Results / Output:**

usr/bin/time outputs:

| SLURM process | Time | %CPU | Memory (kbytes)|
|---|---|---|---|
| Sorting input SAM | .25 | 96% | 934024 |
 

**Issues encountered:**
Need to fix cigar string logic, it is currently not accounting for multiple insertion syntaxes being present within cigar string.
 
**Next steps:**
Fix issue stated above and include log file output code. Once completed, need to test code script. Also need to provide some more verbiage in lab notebook about program and detailed purpose.

**Time spent:**
3.5 hours

### 2026-10-25
 
**Goal for today:**
 Finish deduper script

**Steps performed**
1. Finished CIGAR strand logic and tested with pytest
2. Finished output file code portion
 
**Testing Commands Run for CIGAR string:**
```bash
# unit tests testing decode_cigar() - all passed
assert decode_cigar("72M", 3112307, True) == 3112307
assert decode_cigar("72M", 3112307, False) == 3112379
assert decode_cigar("1S71M", 3123097, True) == 3123096
assert decode_cigar("1S71M", 3123097, False) == 3123168
assert decode_cigar("41M53N31M", 3147363, True) == 3147363
assert decode_cigar("41M53N31M", 3147363, False) == 3147488
assert decode_cigar("68M4S", 3369137, False) == 3369209
assert decode_cigar("5S68M4S", 3369137, False) == 3369209
assert decode_cigar("68M7I4S", 3369137, False) == 3369209
``` 

**Next steps:**
Need to test code on tester SAM file. Maybe consider adding a function to be able to handle both gzipped files and plain .sam file

**Time spent:**
1 hour

### 2026-10-27
 
**Goal for today:**
Test code!!

**Steps performed:**
1. Test code by running on test files and entire sorted SAM file
2. Debug code (step above ended up not going well)
 
**Commands / scripts run:**
```bash

# run on entire sorted SAM file
/usr/bin/time -v pixi run ./balmer_deduper.py \
    -u ../files/STL96.txt \
    -f ../files/sorted_insam.sam \
    -o ../files/sorted_outsam.sam

# run on test file
/usr/bin/time -v pixi run ./balmer_deduper.py \
    -u ../files/STL96.txt \
    -f ../unit_tests/input.sam \
    -o ../unit_tests/test_output.sam
``` 
 
**Results / Output:**
Script works!! 

Output:
```
Summary of duplicate filtering in ../files/sorted_insam.sam file.

Number of header lines: 65

Number of unique reads: 13719048

Number of duplicate reads: 4467362

Number of incorrect UMIs found: 0
```


All my debugging code is in test.py with pytest and I also ran my script on test files within the unit_tests/ directory which ultimately solved my primary issue: bitwise() function to determine strandedness.

**Issues encountered (MANY):**

My code output the following incorrect stats on the first round:
```
Summary of duplicate filtering in ../files/sorted_insam.sam file.

Number of header lines: 65

Number of unique reads: 7667251

Number of duplicate reads: 10519159

Number of incorrect UMIs found: 0
```
All my testing to figure this out took place in test.py with pytest and by running on my test files within the unit_tests/ directory. 
```bash
# key to files in unit_tests/
"""
input.sam: input file for deduper script

output.sam: intended output for comparison with actual output when testing

test_output.sam: output file from testing deduper script

test.py: script for testing functions with pytest
"""
```

Problem #1: My bitwise function was not returning a boolean operator when the strand was reverse (see code below), it was returning a None object which didn't throw an error in my code weirdly enough. It was a simple syntax mistake with my if else statement. This was my main problem.

```python
#### incorrect function

def bitwise(flag:int) -> bool:
    if((flag & 16) != 16):
        forward = True
        return forward


#### corrected function

def bitwise(flag:int) -> bool:
    if((flag & 16) != 16):
        forward = True
    else:
        forward = False
    
    return forward
```
Problem #2: In attempting to solve the problem above, I found that I was adding tuples to my f_strand and r_strand sets incorrectly. I was using the set.update() function when I should have been using the set.add() function. This did change my output numbers a bit, but obviously did not correct it like fixing the bitwise operator problem.

Problem #3: I needed to not convert my chromosome variable to an integer since there are character chromosomes: X, Y, MT

Problem #4: I needed to move my while loop break syntax prior to the line.split() function otherwise I received an error when the program tries to split the last line (which is empty) into various variables.
 
**Time spent:**
2.5 hours


### yyyy-mm-dd
 
**Goal for today:**
 
**Steps performed:**

 
**Commands / scripts run:**
```bash


``` 
 
**Results / Output:**


**Issues encountered:**

 
**Next steps:**
 
**Time spent:**
