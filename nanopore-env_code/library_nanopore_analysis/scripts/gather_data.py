import csv
import pandas as pd
from Bio import SeqIO

# Read input file paths from snakemake
consfile_path = snakemake.input["consfile"]
reffile_path = snakemake.input["reffile"]
cigarfile_path = snakemake.input["cigarfile"]
codebook_path = snakemake.input["inv_barcode_codebook"]

# Read the codebook, cigar data, and parse the fasta files into DataFrames
barcode_codebook = pd.read_csv(codebook_path, delimiter="\t")
cigar_data = pd.read_csv(cigarfile_path, delimiter="\t")

# Function to parse fasta files into a DataFrame
def parse_fasta_to_df(filepath, id_splitter):
    records = []
    for record in SeqIO.parse(filepath, "fasta"):
        record_id = int(record.id.split(id_splitter)[1])
        records.append({"barcodeid": record_id, "sequence": str(record.seq)})
    return pd.DataFrame(records)

# Parse the consensus and reference fasta files
cons_df = parse_fasta_to_df(consfile_path, "_")
ref_df = parse_fasta_to_df(reffile_path, "_")

# Merge the DataFrames on 'barcodeid'
merged_df = barcode_codebook.merge(cigar_data, on="barcodeid")
merged_df = merged_df.merge(cons_df.rename(columns={"sequence": "consensus"}), on="barcodeid")
merged_df = merged_df.merge(ref_df.rename(columns={"sequence": "reference"}), on="barcodeid")

# Reorder columns to match the desired output format
output_df = merged_df[["barcodeid", "barcode", "consensus", "reference", "startpos", "cigar"]].copy()
output_df.rename(columns={"startpos": "alignmentstart"}, inplace=True)

# Write the merged data to the output file
output_df.to_csv(snakemake.output["outputfile"], sep="\t", index=False)



# barcode_codebook = {}
# codebook_data = pd.read_csv(codebook_path, delimiter="\t")
# for _, row in codebook_data.iterrows():
#     barcode_codebook[int(row["barcodeid"])] = row["barcode"]

# cigar_dict = {}
# cigar_data = pd.read_csv(cigarfile_path, delimiter="\t")
# for _, row in cigar_data.iterrows():
#     cigar_dict[int(row["barcodeid"])] = tuple((row["startpos"], row["cigar"]))

# consdict = {
#     int(record.id.split(" ")[0].split("_")[1]): str(record.seq)
#     for record in SeqIO.parse(consfile_path, "fasta")
# }
# refdict = {
#     int(record.id.split("_")[1]): str(record.seq)
#     for record in SeqIO.parse(reffile_path, "fasta")
# }

# output_data = [
#     {
#         "barcodeid": key,
#         "barcode": val,
#         "consensus": consdict.copy()[key],
#         "reference": refdict.copy()[key],
#         "alignmentstart": cigar_dict.copy()[key][0],
#         "cigar": cigar_dict.copy()[key][1],
#     }
#     for key, val in barcode_codebook.items()
# ]

# keys = ["barcodeid", "barcode", "consensus", "reference", "alignmentstart", "cigar"]
# with open(snakemake.output["outputfile"], "w") as outfile:
#     dict_writer = csv.DictWriter(outfile, keys, delimiter="\t")
#     dict_writer.writeheader()
#     dict_writer.writerows(output_data)
