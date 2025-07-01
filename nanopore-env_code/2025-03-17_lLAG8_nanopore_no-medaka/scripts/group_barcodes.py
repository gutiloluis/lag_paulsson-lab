import csv
import ast
import os
import shutil
import numpy as np
import pandas as pd
import re

from Bio import SeqIO

from matplotlib import pyplot as plt


def make_seg_dict(gfafile):
    segment_dict = {}
    with open(gfafile, "r") as infile:
        for line in infile:
            if line[0] == "S":
                splitline = line.split("\t")
                segment_dict[splitline[1]] = splitline[2][:-1]
    return segment_dict


# def generate_reference(cigar, segment_dict):
#     if cigar[0] == "+":
#         traversal = cigar[3].split(">")[1:]
#         ref = "".join(list(map(lambda seg: segment_dict[seg], traversal)))
#     elif cigar[0] == "-":
#         traversal = cigar[3].split("<")[1:][::-1]
#         ref = "".join(list(map(lambda seg: segment_dict[seg], traversal)))
#     return ref


def parse_cigar(cigar):
    
    # Use regex to split the string while keeping delimiters
    components = re.split(r'(>[^<>]+|<[^<>]+)', cigar[3])
    # Filter out any empty strings that may result from the split
    components = [comp for comp in components if comp]
    return components

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(seq))

def generate_reference(cigar, sequence_dict):
    components = parse_cigar(cigar)
    
    reference = ""
    
    for comp in components:
        if comp.startswith('>'):
            key = comp[1:]
            reference += sequence_dict.get(key, comp)
        elif comp.startswith('<'):
            key = comp[1:]
            reversed_seq = reverse_complement(sequence_dict.get(key, comp))
            reference += reversed_seq
        else:
            reference += comp
    if cigar[0] =='-':
        reference = reverse_complement(reference)
    return reference


### Get Barcode to Read Map ###

inv_codebook_path = snakemake.input["inv_barcode_codebook"]
tsv_paths = snakemake.input["tsv"]
fastqs_path = snakemake.input["fastq"]
print(inv_codebook_path)
print(tsv_paths)
print(fastqs_path)
inv_barcode_codebook = {}
codebook_data = pd.read_csv(inv_codebook_path, delimiter="\t")
for _, row in codebook_data.iterrows():
    inv_barcode_codebook[int(row["barcodeid"])] = ast.literal_eval(row["readlist"])

cigar_dict = {}
for filepath in tsv_paths:
    with open(filepath, "r") as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            cigar_dict[row["readname"]] = ast.literal_eval(row["cigar"])

### Make Segment Dictionary ###
gfafile = snakemake.input["gfa"]
seg_dict = make_seg_dict(gfafile)

import psutil

### Make fastq lookup ###
fastqlist = snakemake.input["fastq"]
record_dict = {}
for i, fastqfile in enumerate(fastqlist):
    read_name_list = []
    read_list = []
    with open(fastqfile, "r") as infile:
        for i, line in enumerate(infile):
            if i % 4 == 0:
                read_name_list.append(line.rstrip().split()[0][1:])
            elif i % 4 == 1:
                read_list.append(line.rstrip())
    working_dict = dict(zip(read_name_list, read_list))
    record_dict.update(working_dict)
    del working_dict

readgroup_path = snakemake.output["readgroups"]
if os.path.exists(readgroup_path):
    shutil.rmtree(readgroup_path)
os.makedirs(readgroup_path)

groupref_path = snakemake.output["grouprefs"]
if os.path.exists(groupref_path):
    shutil.rmtree(groupref_path)
os.makedirs(groupref_path)

chunksize = snakemake.params["chunksize"]
subsample_list = snakemake.params["subsample_list"]

for key, val in inv_barcode_codebook.items():

    ref_str = ""
    ref_seq = generate_reference(cigar_dict[val[0]], seg_dict)
    ref_str += ">group_" + str(key) + "\n"
    ref_str += ref_seq + "\n"

    for subsample in subsample_list:
        if not os.path.exists(groupref_path +  "/subsample=" + str(subsample)):
            os.makedirs(groupref_path +  "/subsample=" + str(subsample))
        if not os.path.exists(readgroup_path +  "/subsample=" + str(subsample)):
            os.makedirs(readgroup_path +  "/subsample=" + str(subsample))
        
        
        refchunkpath = (
            groupref_path
            + "/subsample="
            + str(subsample)
            + "/chunk_"
            + str(key // chunksize)
            + ".fasta"
        )
        
        readchunkpath = (
            readgroup_path
            + "/subsample="
            + str(subsample)
            + "/chunk_"
            + str(key // chunksize)
            + ".fastq"
        )

        if subsample == None:
            subsampled_reads = val
        else:
            subsampled_reads = list(
                np.random.choice(val, size=subsample, replace=False)
            )

        out_str = ""
        
        for read_num,read_name in enumerate(subsampled_reads):
            seq_record = record_dict[read_name]
            out_str += "@" + str(key) +"_" +str(read_num)+ "\n"
            out_str += seq_record + "\n"
        
        if key % chunksize == 0 :
            with open(refchunkpath, "w") as outfile:
                outfile.write(ref_str)
            with open(readchunkpath, "w") as outfile:
                outfile.write(out_str)
        
        else:
            with open(refchunkpath, "a") as outfile:
                outfile.write(ref_str)
            with open(readchunkpath, "a") as outfile:
                outfile.write(out_str)
