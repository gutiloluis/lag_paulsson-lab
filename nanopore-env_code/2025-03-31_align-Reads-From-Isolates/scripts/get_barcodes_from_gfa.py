import numpy as np
import pandas as pd

def get_barcode(barcode_list: list) -> str:
    '''
    Takes a list of elements of a cigar path and returns a corresponding 30-mer of 0s and 1s.
    '''
    return "".join(np.array(['ON' in element for element in barcode_list if 'BIT' in element], dtype=int).astype(str))

gaf_columns = [
    "query_name",            # Name of the query sequence
    "query_length",          # Length of the query sequence
    "query_start",           # Start position of the alignment on the query sequence (0-based)
    "query_end",             # End position of the alignment on the query sequence (0-based, exclusive)
    "strand",                # Strand of the alignment ('+' for forward, '-' for reverse)
    "path",                  # Path in the graph that the query aligns to (nodes separated by commas)
    "path_length",           # Length of the path in the graph
    "path_start",            # Start position of the alignment on the path (0-based)
    "path_end",              # End position of the alignment on the path (0-based, exclusive)
    "residue_matches",       # Number of matching bases in the alignment
    "alignment_block_length",# Total length of the alignment block
    "mapping_quality",       # Mapping quality score
    "optional1",               # Optional fields (e.g., tags like NM:i, AS:f, etc.)
    "optional2",              # Optional fields (e.g., tags like NM:i, AS:f, etc.)
    "optional3",              # Optional fields (e.g., tags like NM:i, AS:f, etc.)
    "optional4",              # Optional fields (e.g., tags like NM:i, AS:f, etc.)
    "optional5",              # Optional fields (e.g., tags like NM:i, AS:f, etc.)
]

filename_gaf = snakemake.input[0]
filename_gaf_pickle = snakemake.output[0]

df =(\
    pd.read_csv(filename_gaf,
                sep='\t',
                names = gaf_columns,
                header = None,
                usecols = ['query_name', 'strand', 'path'])#, 'Start', 'End', 'col_15'])
)

df = (df
      .drop_duplicates(subset = 'query_name', keep = False) # This is to discard "chimeric" reads. Maybe they can be handled?
      .astype({'strand': 'category'})#, 'Start': 'int32', 'End': 'int32'}) # Set dtypes
      .loc[lambda df_: (df_['path'].str.contains('SPACER2')) & (df_['path'].str.contains('Variant'))] # Keep only complete reads.
      .assign(
              query_name = lambda df_: df_['query_name'].str.replace(' ', '_'),
              strand = (lambda df_: df_['strand']
                        .cat.set_categories(['+', '-'])
                        .mask(df_['path'].str.contains('<'), '-')
                       )
             )
     )

# TODO: Figure out how to chain these two parts
df.loc[df['strand'] == '+', 'barcode'] = (df
                                          .loc[df['strand'] == '+','path']
                                          .str.split('>')
                                          .apply(get_barcode)
                                         )# Get barcode list

df.loc[df['strand'] == '-', 'barcode'] = (df.
                                          loc[df['strand'] == '-','path']
                                          .str.split('<').str[::-1]
                                          .apply(get_barcode)
                                         )# Get barcode list

df.to_pickle(filename_gaf_pickle) # Save to pickle