import ast
import os
import shutil
import numpy as np
import pandas as pd

## Get Segments to sequences from graph reference
with open(snakemake.input['gfa'], 'r') as infile:
    x = np.array(infile.readlines())
    
last_row_of_seq = np.argwhere(np.diff([not line.startswith('S') for line in x[1:]]))[0][0] + 1

df_segs_to_seq = (pd
                  .read_csv(snakemake.input['gfa'],
                            sep='\t',
                            nrows= last_row_of_seq,
                            skiprows = 1,
                            header=None,
                            names = [0, 'seg_name', 'seg_sequence'],
                            usecols=['seg_name', 'seg_sequence'],
                            index_col='seg_name')
                  .iloc[:,0]
                 )

#######################
## Read all the tsvs with the read ids and paths. This is faster than loading the cigar dict!
inpathlist = snakemake.input["tsv"]
# Read all dataframes and concatenate them
df_reads = pd.concat((pd
                      .read_csv(filename,
                                sep = '\t',
                                usecols=['Query_name', 'Barcode', 'Strand', 'Path_matching'])
                      .astype({'Strand':'category'})
                     
                for filename in inpathlist),
               ignore_index = True)

########################
## Get barcode df. Each cell in the 'Query_name' col contains the list of reads corresponding to a given bc
df_codebook = pd.read_csv(snakemake.input['inv_barcode_codebook'], sep = '\t')

df_codebook = (df_codebook
               .astype({'Barcode_id': 'int32', 'Counts': 'int32'})
               .assign(Query_name = lambda df_:df_['Query_name']
                       .astype(str)
                       .apply(ast.literal_eval))
               .assign(First_read = lambda df_:df_['Query_name'].str[0])
              )

# Merge read seq
df_codebook = (df_codebook
               .merge(df_reads[['Query_name', 'Strand', 'Path_matching']],
                      left_on='First_read',
                      right_on='Query_name',
                      how='left',
                      suffixes=('','_y'),
                      validate='one_to_many')
               .drop(columns=['Query_name_y'])
              )

# Get seg_list # TODO: DO IT ON THE FIRST get_barcodes.py. FASTEST WAY?
df_codebook.loc[df_codebook['Strand'] == '+', 'Seg_list'] = df_codebook.loc[df_codebook['Strand'] == '+', 'Path_matching'].str.split('>').str[1:]
df_codebook.loc[df_codebook['Strand'] == '-', 'Seg_list'] = df_codebook.loc[df_codebook['Strand'] == '-', 'Path_matching'].str.split('<').str[1:].str[::-1]

def generate_reference(seg_list): # Not pass df_segs_to_seq
    return ''.join(list(map(lambda segment: df_segs_to_seq.loc[segment], seg_list)))

# Get ref_seq
df_codebook['ref_seq'] = df_codebook['Seg_list'].apply(generate_reference)


#############
# df having name of read and read sequence (from fastq)
def parse_fastq_files(fastqlist):
    '''
    fastqlist (list): List of fastq filenames with the reads
    
    Returns
        pd.DataFrame: df with read_name as index and read_sequence
    '''
    record_dict = {}
    for i, fastqfile in enumerate(fastqlist):
        read_name_list = []
        read_list = []
        with open(fastqfile, "r") as infile:
            for i, line in enumerate(infile):
                if i % 4 == 0:
                    read_name_list.append(
                        line.rstrip().replace(" ", "_")[1:]
                    )  ##changed to include whole readname, whitespace replaced with underscore
                elif i % 4 == 1:
                    read_list.append(line.rstrip())
        working_dict = dict(zip(read_name_list, read_list))
        record_dict.update(working_dict)
        del working_dict
    
    return pd.DataFrame.from_dict(record_dict, orient='index', columns=['Read_sequence'])

df_reads = parse_fastq_files(snakemake.input['fastq'])

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


# TODO INCLUDE SUBSAMPLING, PASS DFS
def write_reference_to_file(barcode_id, ref_seq, read_id_list):
    
    # Get string for reference (fasta)
    ref_str = ''
    ref_str = '>group_' + str(barcode_id) + '\n'
    ref_str += ref_seq + '\n'

    subsample = None

    intermediate_path = (
        '/subsample='
        + str(subsample)
        + '/chunk_'
        + str(barcode_id // chunksize)   
    )
    
    refchunkpath = (
        groupref_path
        + intermediate_path
    )

    if not os.path.exists(refchunkpath):
                os.makedirs(refchunkpath)

    ref_fasta_path = refchunkpath + '/group_' + str(barcode_id) + '.fasta'

    with open(ref_fasta_path, 'w') as outfile:
        outfile.write(ref_str)
        
    readchunkpath = (
        readgroup_path
        + intermediate_path
    )
    
    # if else subsampel
    if not os.path.exists(readchunkpath):
        os.makedirs(readchunkpath)
    
    read_str = ''
    for read_id in read_id_list:
        read_seq = df_reads.loc[read_id, 'Read_sequence']
        read_str += '@' + read_id + '\n'
        read_str += read_seq + '\n'
    
    with open(readchunkpath + '/group_' + str(barcode_id) + '.fastq', 'w') as outfile:
        outfile.write(read_str)

# Write files for each barcode (a reference file and all the reads)
_ = df_codebook.apply(lambda df_: write_reference_to_file(df_['Barcode_id'], df_['ref_seq'], df_['Query_name']), axis=1)