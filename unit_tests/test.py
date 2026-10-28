#!/usr/bin/env python

import re

def decode_cigar(cigar: str, raw_stt: int, str: bool) -> int:
    '''takes in cigar string, start site, and strand (T = forward, F = reverse) to returns true start site accounting for soft clipping on 5' end'''
    numbers = re.findall("[0-9]+", cigar)
    letters = re.findall("[MNIDS]", cigar)

    if str:
        if letters[0] == "S":
            true_stt = raw_stt - int(numbers[0])
        else:
            true_stt = raw_stt

    else:
        for pos, l in enumerate(letters):

            if pos == 0 and l == "S":
                continue
            
            elif l == "I":
                continue

            raw_stt += int(numbers[pos])

        true_stt = raw_stt
        
        
    return true_stt


f_strand = set()
r_strand = set()

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

def bitwise(flag:int) -> bool:
    '''returns strandedness via Boolean value based on bitwise flag (True = forward, False = reverse)'''
    
    if((flag & 16) != 16):
        forward = True
    else:
        forward = False
    
    return forward
    
'''
print(uniq_check(76718921, True, "ATCGAACC"))
print(f_strand)
print(r_strand)

print(uniq_check(76718921, True, "ATCGAACC"))
print(f_strand)
print(r_strand)

print(uniq_check(76718921, False, "ATCGAACC"))
print(f_strand)
print(r_strand)

print(uniq_check(76718921, False, "ATCGAACC"))
print(f_strand)
print(r_strand)
'''

def test_bitwise():
    assert bitwise(16) == False
    assert bitwise(0) == True


def test_decode_cigar():
    assert decode_cigar("74M", 76718921, True) == 76718921
    assert decode_cigar("3S71M", 76718924, True) == 76718921
    assert decode_cigar("1S71M", 3123097, True) == 3123096
    assert decode_cigar("1S71M", 3123097, False) == 3123168
    assert decode_cigar("41M53N31M", 3147363, True) == 3147363
    assert decode_cigar("41M53N31M", 3147363, False) == 3147488
    assert decode_cigar("68M4S", 3369137, False) == 3369209
    assert decode_cigar("5S68M4S", 3369137, False) == 3369209
    assert decode_cigar("68M7I4S", 3369137, False) == 3369209
    assert decode_cigar("4S25M470N46M", 76875953, True) == 76875949
    assert decode_cigar("6S15M470N46M", 76875955, True) == 76875949



def test_uniq_check():
    assert uniq_check(76718921, True, "ATCGAACC") == True
    assert uniq_check(76718921, True, "ATCGAACC") == False
    assert uniq_check(76718921, False, "ATCGAACC") == True
    assert uniq_check(76718921, False, "ATCGAACC") == False




