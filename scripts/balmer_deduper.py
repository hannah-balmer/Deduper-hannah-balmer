#!/usr/bin/env python

import re
import argparse

def get_args():
    parser = argparse.ArgumentParser(description = "Identify and filter out PCR duplicate reads from a single-end sequencing post-alignment SAM file. Will output a filtered SAM file and a .log file containing metrics on duplicate removal.")
    parser.add_argument("-u", "--umi", help="Absolute file path to a text file containing UMI sequences newline separated.", type = str)
    parser.add_argument("-f", "--file", help="Absolute path of the single-end sequencing post-alignment SAM file. Must be sorted by chromsome then start site.", type = str)
    parser.add_argument("-o", "--outfile", help="Absolute file path where you would like the duplicate-filtered SAM file to be created.", type = str)
    return parser.parse_args()

args = get_args()

txt_umi = args.umi
in_sam = args.file
out_sam = args.outfile


def decode_cigar(cigar: str, raw_stt: int, str: bool) -> int:
    '''takes in cigar string, start site, and strand (T = forward, F = reverse) to returns true start site accounting for soft clipping on 5' end'''
    
    # regex to isolate all letters and numbers within CIGAR string to separate lists
    numbers = re.findall("[0-9]+", cigar)
    letters = re.findall("[MNIDS]", cigar)

    # if forward strand, adjust for 5' soft clipping is present
    if str:
        if letters[0] == "S":
            true_stt = raw_stt - int(numbers[0])
        else:
            true_stt = raw_stt

    # if reverse strand, add all cigar string numbers to start site except for 3' soft clipping and any insertions
    else:
        for pos, l in enumerate(letters):

            if pos == 0 and l == "S":
                continue
            
            elif l == "I":
                continue

            raw_stt += int(numbers[pos])

        true_stt = raw_stt
        
    return true_stt


def bitwise(flag:int) -> bool:
    '''returns strandedness via Boolean value based on bitwise flag (True = forward, False = reverse)'''
    
    if((flag & 16) != 16):
        forward = True
    else:
        forward = False
    
    return forward

def uniq_check(adj_stt: int, str: bool, umi: str) -> bool:
    '''takes in adjusted start site (output from decode_cigar), strand, and umi to output Boolean value (T for unique read, F for duplicate read).'''
    
    # if read is in forward strand and is not a duplicate
    if str:
        if (adj_stt, umi) not in f_strand:
            f_strand.add((adj_stt, umi))
            return True
        
        # if read is a duplicate
        else:
            return False

    # if read is in reverse strand and is not a duplicate
    elif str == False:
        if (adj_stt, umi) not in r_strand:
            r_strand.add((adj_stt, umi))
            return True
        
        # if read is a duplicate
        else:
            return False



# read in UMIs into set
umis = set()
with open(txt_umi, "r") as umifile:
    while True:
        line = umifile.readline().strip()

        if line == "":
            break

        umis.add(line)



# initialize all variables
n_headers = 0 # number of header lines
unique = 0 # number of unique reads
dupes = 0 # number of duplicates filtered out
xUMIs = 0 # number of reads with incorrect UMIs

current_chrom = "1" # current chromosome
f_strand = set() # will hold all UMI and start site combos in tuple for each unique read found in forward strand 
r_strand = set() # will hold all UMI and start site combos in tuple for each unique read found in reverse strand


with open(out_sam, "w") as osam:
    
    with open(in_sam, "r") as isam:
        
        while True:
            line = isam.readline().strip()

            # if header line: write to output file 
            if line.startswith("@"):
                n_headers += 1
                osam.write(f"{line}\n")
                continue

            # end loop when end of file reached
            elif line == "":
                break

            # if read line
            else:
                #nline += 1
                line_cols = line.split()

                # isolate necessary columns in read
                qname = line_cols[0]
                bit = int(line_cols[1])
                chrom = line_cols[2]
                start = int(line_cols[3])
                cigar = line_cols[5]

                # isolate UMI
                match = re.search(r":([ATCG]+)$", qname)
                umi = match.group(1)

                # continue to next line if UMI does not match
                if umi not in umis:
                    xUMIs += 1
                    continue

                # to save memeory, clear both strand sets and set chromosome variable to new number value
                if chrom != current_chrom:
                    f_strand.clear()
                    r_strand.clear()
                    current_chrom = chrom

                # determine strand
                strand = bitwise(bit)

                # determine true start site
                Tstart = decode_cigar(cigar, start, strand)

                # determine whether unique or duplicate
                if uniq_check(Tstart, strand, umi):
                    unique += 1
                    osam.write(f"{line}\n")

                else:
                    dupes += 1
                    continue

# write to slurm script 
print(f"Summary of duplicate filtering in {in_sam} file.\n")
print(f"Number of header lines: {n_headers}\n")
print(f"Number of unique reads: {unique}\n")
print(f"Number of duplicate reads: {dupes}\n")
print(f"Number of incorrect UMIs found: {xUMIs}\n")