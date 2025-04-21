import pandas as pd

filename_grnas = snakemake.input[0]
filename_nanopore = snakemake.input[1]
filename_summary = snakemake.output[0]
filename_error_report = snakemake.output[1]

df_isolates = pd.read_pickle(filename_grnas)
df_nanopore = pd.read_pickle(filename_nanopore)

#  Map isolates to nanopore data
df_summary = \
(df_isolates
 [['query_name', 'barcode', 'isolate_number', 'grna_sequence_isolates']]
 .sort_values(by='isolate_number')
 .merge(
    df_nanopore[['barcode', 'closest_match', 'sgRNA Sequence', 'barcode depth']].rename(columns = {'sgRNA Sequence': 'grna_sequence_nanopore', 'closest_match': 'closest_hamming_distance'}),
    on='barcode',
    # Set the suffixes
    how='left'
 )
 .assign(is_barcode_in_nanopore = lambda df_: df_['barcode'].isin(df_nanopore['barcode']))
 .assign(is_grna_sequence_match = lambda df_: df_['grna_sequence_isolates'] == df_['grna_sequence_nanopore'])
)

df_summary.to_pickle(filename_summary)

# Generate summary of results
number_of_isolates = len(df_summary)

summary_lines = []
summary_lines.append('Filename Nanopore: {}'.format(filename_nanopore))
summary_lines.append('Filename grnas: {}'.format(filename_grnas))
summary_lines.append('Number of barcodes in nanopore data: {}'.format(len(df_nanopore)))
summary_lines.append('Minimum barcode depth: {}'.format(df_nanopore['barcode depth'].min()))
summary_lines.append('Number of isolates: {}'.format(number_of_isolates))
number_of_mapped_isolates = df_summary['is_barcode_in_nanopore'].sum()
fraction_of_mapped_isolates = number_of_mapped_isolates / number_of_isolates
summary_lines.append('Number of isolates mapped to nanopore data: {} ({:.2f} %)'.format(number_of_mapped_isolates,
                                                               fraction_of_mapped_isolates*100))

number_of_isolate_barcodes_in_nanopore = len(df_summary[(df_summary['is_barcode_in_nanopore'] == True)])
number_of_mapped_matching_isolates = len(df_summary[(df_summary['is_barcode_in_nanopore'] == True) & (df_summary['is_grna_sequence_match'] == True)]) 
fraction_of_mapped_matching_isolates = number_of_mapped_matching_isolates / number_of_isolate_barcodes_in_nanopore if number_of_isolate_barcodes_in_nanopore != 0 else 0
summary_lines.append('Number of isolates mapped to nanopore data with matching gRNA sequences: {} ({:.2f} %)'.format(number_of_mapped_matching_isolates,
                                                                                    fraction_of_mapped_matching_isolates*100))
summary_lines.append('\n')

# Print the summary
for line in summary_lines:
   print(line)

# Save the summary to the error report file
with open(filename_error_report, 'w') as f:
   f.write('\n'.join(summary_lines))