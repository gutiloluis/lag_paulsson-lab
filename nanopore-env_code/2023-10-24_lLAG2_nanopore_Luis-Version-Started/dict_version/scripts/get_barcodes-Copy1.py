import csv
import numpy as np

cigar_dict = {}

snakemake_input = '/home/lag36/scratch/lag36/2023-10-24_lDE2-nanopore/graph_output/tube_1_01.gaf'
snakemake_output = '/home/lag36/scratch/lag36/2023-10-24_lDE2-nanopore/graph_output/tube_1_01.tsv'

with open(snakemake_input, "r") as infile:
    for line in infile:
        data = line.split("\t")
        ##changed to include whole readname, whitespace replaced with underscore
        read_id = data[0].replace(" ", "_")
        if ">" in data[5]:
            cigar_dict[read_id] = (
                "+",
                int(data[7]),
                int(data[8]),
                data[5],
                data[15].split(":")[-1][:-1],
            )
        else:
            cigar_dict[read_id] = (
                "-",
                int(data[7]),
                int(data[8]),
                data[5],
                data[15].split(":")[-1][:-1],
            )

barcode_dict = {}
for key in cigar_dict.keys():
    cigar = cigar_dict[key]
    if "Variant" in cigar[3] and "SPACER2" in cigar[3]:
        if cigar[0] == "-":
            barcode = cigar[3].split("<")
            barcode = barcode[::-1]
        elif cigar[0] == "+":
            barcode = cigar[3].split(">")
        barcode = barcode[:-1]
        barcode = (
            np.array(["ON" in item for item in barcode if "BIT" in item])
            .astype(int)
            .astype(str)
            .tolist()
        )
        barcode = "".join(barcode)
        barcode_dict[key] = barcode

barcode_dict = [
    {"readname": key, "barcode": val, "cigar": cigar_dict[key]}
    for key, val in barcode_dict.items()
]

keys = ["readname", "barcode", "cigar"]

# print(barcode_dict)

with open(snakemake_output, "w") as outfile:
    dict_writer = csv.DictWriter(outfile, keys, delimiter="\t")
    dict_writer.writeheader()
    dict_writer.writerows(barcode_dict)
