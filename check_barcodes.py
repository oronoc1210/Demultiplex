import subprocess
import os
import re

import pandas as pd

def barcode_count(library, barcode):
    cat_cmd = subprocess.Popen(('cat', library), stdout=subprocess.PIPE)
    grep_cmd = subprocess.Popen(('grep', barcode), stdin=cat_cmd.stdout, stdout=subprocess.PIPE)
    output = subprocess.check_output(('wc', '-l'), stdin=grep_cmd.stdout)
    return output

def main():
    counting_dict = dict()
    barcode_sequences = {'CGGACTTCTGTA': 13,
                         'CATATGGAACCG': 14,
                         'GCACACCTATAC': 15,
                         'ACACTTGGCCTC': 16,
                         'TTCATAACGCCA': 17,
                         'GCACTCGTAACT': 18,
                         'TTCGATCAATCC': 19,
                         'ACGGCGAGTTAT': 20,
                         'GGTTGAGGATCA': 21,
                         'TCATATGGCGCG': 22,
                         'GATCTGCGGTGT': 23,
                         'AGATCCTTAGAG': 24}
    counting_dict['Barcode Sequence'] = dict()
    for barcode_sequence in barcode_sequences.keys():
        barcode_id = "Barcode {}".format(barcode_sequences[barcode_sequence])
        counting_dict['Barcode Sequence'][barcode_id] = barcode_sequence
    for file in os.listdir(os.getcwd()):
        if re.match(r'([A-Z]{5}).fastq', file) is not None:
            print("Parsing library {}...".format(file))
            counting_dict[file] = dict()
            for barcode_sequence in barcode_sequences.keys():
                count = barcode_count(file, barcode_sequence)
                barcode_id = "Barcode {}".format(barcode_sequences[barcode_sequence])
                counting_dict[file][barcode_id] = count
    print("Creating excel table...")
    df = pd.DataFrame.from_dict(counting_dict)
    print(df)
    df.to_excel("barcode_estimates.xlsx")




if __name__ == '__main__':
    main()
