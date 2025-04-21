from Bio import SeqIO
import pandas as pd
from Bio.Seq import Seq

PROTOSPACER_UPSTREAM_SEQUENCE = 'TTTCTAGCTCTAAAAC' # If I apply it to other backbone designs, I should change this
HEADPATH = '/home/lag36/scratch/lag36/2025-03-31_align-Reads-From-Isolates/'

# inputs
filename_fasta = snakemake.input[0]
filename_gaf_pickle = snakemake.input[1]
filename_grnas = snakemake.output[0]

def read_fasta(filename):
    fasta_reads = {}
    for fasta_record in SeqIO.parse(filename, "fasta"):
        read_number = int(fasta_record.name
                          .split('_')[2]
                          .split('-')[1])
        read_sequence = fasta_record.seq
        fasta_reads[read_number] = read_sequence
    # Convert to DataFrame
    fasta_reads = \
    (pd.DataFrame
     .from_dict(fasta_reads, orient='index', columns=['sequence'])
     .reset_index()
     .rename(columns={'index': 'isolate_number'})
    )
    return fasta_reads

# Read fasta file
reads_isolates = read_fasta(filename_fasta)

# Read and process the isolates data
df_isolates = \
(pd.read_pickle(filename_gaf_pickle)
 
 # Parse isolate number from filename. THIS NEEDS TO BE MODIFIED IF FILENAME FORMAT CHANGES
 .assign(isolate_number=lambda df_: df_['query_name']
         .str.split('_')
         .str[2]
         .str.split('-')
         .str[1]
         .astype(int))
 .merge(reads_isolates,
        on='isolate_number',
        how='left')
 .assign(sequence = lambda df_: df_.apply(lambda row:
                                          row['sequence'].reverse_complement() if row['strand']=='-' else row['sequence'], axis=1
                                          ).astype(str))
 .assign(protospacer_start = lambda df_: df_['sequence'].str.find(PROTOSPACER_UPSTREAM_SEQUENCE) + len(PROTOSPACER_UPSTREAM_SEQUENCE))
 .assign(protospacer_sequence = lambda df_: df_.apply(lambda row:
                                             row['sequence'][row['protospacer_start']:row['protospacer_start']+20], axis=1
                                             ))
 .assign(grna_sequence_isolates = lambda df_: df_['protospacer_sequence'].apply(lambda seq: str(Seq(seq).reverse_complement())))
 .drop(columns=['protospacer_sequence'])
)

df_isolates.to_pickle(filename_grnas)