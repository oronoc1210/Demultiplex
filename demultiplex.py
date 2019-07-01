import os
import sys
import re
import subprocess
import time
import collections
import contextlib

import HTSeq

def test_lib_to_regex(lib2bar, bar2reg):
    """Tests to make sure the correct regex sequences 
    are selected for each library"""
    for lib in lib2bar.keys():
        sequences = [seq for seq in bar2reg.keys() if any(num in seq for num in lib2bar[lib])]
        #sequences = [seq if any(num in seq for num in lib2bar[lib]) else pass for seq in bar2reg.keys()]
        print("{}: {}".format(lib, sequences))

def get_counts(sequences, infile, seq2regex):
    """Get number of reads within a file that contain
    each sequence (sequences should be a list)"""
    fastq_file = HTSeq.FastqReader(infile)
    counts = collections.Counter()
    old_id = None
    old_seq = None
    barcode = False
    for read in fastq_file:
        counts["reads"] += 1
        read_id = read.name.split(' ')[0]
        if old_id == read_id:
            counts["read pairs"] += 1
        else:
            old_id = read_id
            if barcode == False:
                counts["nothin"] += 1
            barcode = False
            old_seq = read.seq
        for sequence in sequences:
            if re.search(seq2regex[sequence], read.seq):
                if barcode:
                    counts["double"] += 1
                if barcode==sequence:
                    counts["double same"] += 1
                counts[sequence] += 1
                barcode=sequence
    print("Reads: %d" % (counts["reads"]))
    print("Read pairs: %d" % (counts["read pairs"]))
    print("Nothin: %d" % (counts["nothin"]))
    for sequence in sequences:
        print("%s: %d" % (sequence, counts[sequence]))
    print("Double: %d" % (counts["double"]))
    print("Double same: %d" % (counts["double same"]))
    print("\n")

def demultiplex(infile, outfile, sequences, seq2regex):
    fastq_file = HTSeq.FastqReader(infile)
    with open(outfile, "w+") as outf:
        for read in fastq_file:
            for sequence in sequences:
                if 'r' not in sequence:
                    continue
                match = re.search(seq2regex[sequence], read.seq)
                if not match:
                    continue
                barcode = HTSeq.Sequence(match.group(0))
                read2 = read.trim_left_end_with_quals(barcode)
                read2.write_to_fastq_file(outf)
                break


def separate_demultiplex(libraries_dir, libname, libpath, lib2seq, sequences, seq2regex):
    fastq_file = HTSeq.FastqReader(libpath)
    seq_storage = collections.Counter()
    count = 0
    for barcode_num in lib2seq[libname]:
        seq_storage[barcode_num] = []
    for read in fastq_file:
        for sequence in sequences:
            if 'r' not in sequence:
                continue
            match = re.search(seq2regex[sequence], read.seq)
            if not match:
                continue
            count += 1
            barcode = HTSeq.Sequence(match.group(0))
            barcode_num = sequence[:2]
            read2 = read.trim_left_end_with_quals(barcode)
            seq_storage[barcode_num].append(read2)
            break
    for barcode_num in seq_storage.keys():
        with open("{}/{}/{}.filter-RNA.demulti.{}.fastq".format(libraries_dir, libname, libname, barcode_num), "w+") as outf:
            for read in seq_storage[barcode_num]:
                read.write_to_fastq_file(outf)

def main():
    lib2seq = {"GGCWY": ["13"],
               "GGCWZ": ["14"],
               "GGCXA": ["15"],
               "GGCXB": ["16"],
               "GGCXC": ["17"],
               "GGCXG": ["18"],
               "GGCXH": ["19"],
               "GGCXN": ["17", "18", "19"],
               "GGCXO": ["13", "14", "15", "16"]
               }

    seq2regex = {"13f": 'CGGACTTCTGTA.*',
                 "13r": '.*TACAGAAGTCCG',
                 "14f": 'CATATGGAACCG.*',
                 "14r": '.*CGGTTCCATATG',
                 "15f": 'GCACACCTATAC.*',
                 "15r": '.*GTATAGGTGTGC',
                 "16f": 'ACACTTGGCCTC.*',
                 "16r": '.*GAGGCCAAGTGT',
                 "17f": 'TTCATAACGCCA.*',
                 "17r": '.*TGGCGTTATGAA',
                 "18f": 'GCACTCGTAACT.*',
                 "18r": '.*AGTTACGAGTGC',
                 "19f": 'TTCGATCAATCC.*',
                 "19r": '.*GGATTGATCGAA'
                }
                     
    libraries_dir = "/global/projectb/scratch/cmodonog/Multiplex_barcoding/Libraries"
    lib1 = "{}/GGCXO/GGCXO.filter-RNA.fastq".format(libraries_dir)
    lib2 = "{}/GGCXN/GGCXN.filter-RNA.fastq".format(libraries_dir)
    libs = [lib1, lib2]
    #for lib in lib2seq.keys():
    for lib in libs:
        libname = re.search('([A-Z]{5})', lib).group(1)
        sequences = [seq for seq in seq2regex.keys() if any(num in seq for num in lib2seq[libname])]
        print("Demultiplexing library {}".format(lib))
        separate_demultiplex(libraries_dir, libname, lib, lib2seq, sequences, seq2regex)
        """
        infile = "{}/{}/{}.filter-RNA.fastq".format(libraries_dir, lib, lib)
        outfile = "{}/{}/{}.filter-RNA.demulti2.fastq".format(libraries_dir, lib, lib)
        print("Demultiplexing library {}".format(lib))
        demultiplex(infile, outfile, sequences, seq2regex)
        """
        #infile = "{}/{}/{}.filter-RNA.fastq".format(libraries_dir, lib, lib)
        print(lib)
        get_counts(sequences, lib, seq2regex)

if __name__ == "__main__":
    main()
