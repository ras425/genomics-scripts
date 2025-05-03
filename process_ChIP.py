#!/usr/bin/env python3

'''
Example usage:
chmod +x process_ChIP.py

./process_ChIP.py text-to-bed WT_vs_KO_narrowPeak.txt WT_vs_KO.bed
./process_ChIP.py extend-reads WT_vs_KO.bed
./process_ChIP.py binding-sites WT_vs_KO.bed WT_vs_KO_sites.bed

Paste into GREAT ChIP, view all region-gene associations, download gene -> genomic region association table
./process_ChIP.py process WT_vs_KO-all-gene.txt WT_vs_KO-2k-genes.txt

'''

import csv
import re
import pandas as pd 
import argparse 
import tempfile
import os

def parse_args():
    parser = argparse.ArgumentParser(description="ChIP-seq processing pipeline")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # text-to-bed
    parser_text = subparsers.add_parser("text-to-bed", help="Convert NarrowPeak text file to BED format")
    parser_text.add_argument("input", help="Input text file")
    parser_text.add_argument("output", help="Output BED file")

    # extend-reads
    parser_extend = subparsers.add_parser("extend-reads", help="Extend reads in a BED file shorter than a specified threshold")
    parser_extend.add_argument("input", help="Input BED file")
    parser_extend.add_argument("output", help="Output BED file with extended reads", nargs="?", default=None)
    parser_extend.add_argument("--threshold", type=int, default=200, help="Length threshold to extend reads (default: 200)")

    # binding-sites
    parser_binding = subparsers.add_parser("binding-sites", help="Extract top N binding sites from BED file")
    parser_binding.add_argument("input", help="Input BED file")
    parser_binding.add_argument("output", help="Output BED file with top binding sites")
    parser_binding.add_argument("--lines", type=int, default=10000, help="Number of lines to extract (default: 10000)")

    # process gene list
    parser_process = subparsers.add_parser("process", help="Extract top genes from GREAT ChIP output")
    parser_process.add_argument("input", help="Input GREAT ChIP output file")
    parser_process.add_argument("output", help="Output file for top genes")
    parser_process.add_argument("--top-n", type=int, default=2000, help="Number of top genes to extract (default: 2000)")

    return parser.parse_args()

def text_to_bed(textfile, bedout):
    '''
    Converts a NarrowPeak format text file to BED format
    '''
    try:
        with open(textfile, "r") as infile, open(bedout, "w", newline="") as outfile:
            writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
            
            for line_num, line in enumerate(infile):
                if line_num == 0:
                    continue  # skip header line
                
                row = line.strip().split()
                chromosome = f"chr{row[0]}"
                start, end = row[1], row[2]

                writer.writerow([chromosome, start, end])
    except Exception as e:
        print(f"Error in text_to_bed: {e}")

def extend_reads(bedfile, output_file, threshold=200):
    '''
    Extends reads from a BED file that are shorter than the specified threshold
    If no output file specified, overwrites the input file!
    '''
    try:
        input_handle = open(bedfile, "r")
        if output_file:
            output_handle = open(output_file, "w", newline="")
            tmp_path = None
        else:
            tmp = tempfile.NamedTemporaryFile(mode='w', delete=False, newline="")
            output_handle = tmp
            tmp_path = tmp.name

        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")

        for line in input_handle:
            row = line.strip().split()

            chromosome = row[0]
            start, end = int(row[1]), int(row[2])

            read_length = end - start
            if read_length < threshold:
                midpoint = (start + end) // 2
                new_start = max(midpoint - threshold // 2, 0)
                new_end = new_start + threshold
                writer.writerow([chromosome, new_start, new_end])
            else:
                writer.writerow([chromosome, start, end])

        input_handle.close()
        output_handle.close()

        if tmp_path:
            os.replace(tmp_path, bedfile)

    except Exception as e:
        print(f"Error in extend_reads: {e}")


def binding_sites(bedfile, output, lines=10000):
    '''
    Extract the top 'lines' number of binding sites from a BED file
    '''
    try:
        df = pd.read_csv(bedfile, sep='\t', header=None, usecols=[0, 1, 2])
        df.head(lines).to_csv(output, header=False, index=False, sep='\t')
    except Exception as e:
        print(f"Error in binding_sites: {e}")

def process(filename, output, top_n=2000):
    '''
    Processes Great ChIP output to get the top 'top_n' genes closest to the TSS
    '''
    try:
        with open(filename, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            next(reader)  # skip header
            
            processed_rows = []
            for row in reader:
                if row:
                    match = re.search(r'\(([-+]?\d+)\)', row[1])
                    if match:
                        number = abs(int(match.group(1)))
                        processed_rows.append([row[0], number])
            
            # sort rows by distance from TSS
            processed_rows.sort(key=lambda x: x[1])
            
            with open(output, 'w') as output_file:
                for i, processed_row in enumerate(processed_rows[:top_n]):
                    output_file.write(f"{processed_row[0]}\n")
    except Exception as e:
        print(f"Error in process: {e}")

def main():
    args = parse_args()

    if args.command == "text-to-bed":
        text_to_bed(args.input, args.output)
    elif args.command == "extend-reads":
        extend_reads(args.input, args.output, args.threshold)
    elif args.command == "binding-sites":
        binding_sites(args.input, args.output, args.lines)
    elif args.command == "process":
        process(args.input, args.output, args.top_n)

if __name__ == "__main__":
    main()