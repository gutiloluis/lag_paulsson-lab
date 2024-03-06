# TODO: What is 15?
# TODO: Use chaining to make it cleaner
# TODO: Use the appropriate data types
import numpy as np
import pandas as pd

def get_barcode(barcode_list):
    '''
    Takes a list of elements of Path_matching and returns a corresponding 30-mer of 0s and 1s.
    '''
    return "".join(np.array(['ON' in element for element in barcode_list if 'BIT' in element], dtype=int).astype(str))

col_names = ['Query_name',
          'Query_length',
          'Query_start',
          'Query_end',
          'Strand',
          'Path_matching',
          'Path_length',
          'Start',
          'End',
          'Number_of_residue_matches',
          'Alignment_block_length',
          'Mapping_quality',
          'col_12',
          'col_13',
          'col_14',
          'col_15',
          'col_16']

snakemake_input = '/home/lag36/scratch/lag36/2023-10-24_lDE2-nanopore/graph_output/tube_1_01.gaf'
snakemake_output = '/home/lag36/scratch/lag36/2023-10-24_lDE2-nanopore/graph_output/tube_1_01_pd.tsv'

df = pd.read_csv(snakemake.input[0],
                  sep='\t',
                  names = col_names,
                  header = None,
                  usecols = ['Query_name', 'Strand', 'Path_matching'])#, 'Start', 'End', 'col_15'])
     
df = (df
      .drop_duplicates(subset = 'Query_name', keep = False) # This is to discard "chimeric" reads. Maybe they can be handled?
      
      .astype({'Strand': 'category'})#, 'Start': 'int32', 'End': 'int32'}) # Set dtypes
      
      .query("Path_matching.str.contains('SPACER2') and Path_matching.str.contains('Variant')") # Keep only complete reads. TODO Could this cause problems on future calls when chaining? MAYBE I NEED TO USE FUNCTIONS!
      
      .assign(#col_15 = df['col_15'].astype(str).str.split(':').str[-1].astype('float32'), # Format col 15
              Query_name = lambda df_: df_['Query_name'].str.replace(' ', '_'),
              # Query_name = df['Query_name'].astype(str).str.replace(' ', '_'),
              Strand = (lambda df_: df_['Strand']
                        .cat.set_categories(['+', '-'])
                        .mask(df_['Path_matching'].str.contains('<'), '-')
                       )
             )
     )

# Get barcodes. TODO: How to incorporate the logic on chain? 
# TODO CONSIDER SAVING THE LIST (path split) HERE! I AM REPEATING THE COMPUTATION TWICE
df.loc[df['Strand'] == '+', 'Barcode'] = (df
                                          .loc[df['Strand'] == '+','Path_matching']
                                          .str.split('>')
                                          .apply(get_barcode)
                                         )# Get barcode list

df.loc[df['Strand'] == '-', 'Barcode'] = (df.
                                          loc[df['Strand'] == '-','Path_matching']
                                          .str.split('<').str[::-1]
                                          .apply(get_barcode)
                                         )# Get barcode list

# df = df.loc[:, ['Query_name', 'Barcode', 'Strand', 'Start', 'End', 'Path_matching', 'col_15']] # Rearrange to match order of dictionary
# df = df.loc[:, ['Query_name', 'Barcode', 'Strand', 'Path_matching']] # Rearrange to match order of dictionary

df.to_csv(snakemake.output[0],
          sep='\t',
          index = False)