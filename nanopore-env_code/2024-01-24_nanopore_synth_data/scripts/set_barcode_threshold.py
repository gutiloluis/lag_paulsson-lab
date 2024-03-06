import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

inpathlist = snakemake.input["tsv"]
# Read all dataframes and concatenate them
# TODO Possible memory issue here if too many files are opened
df = pd.concat((pd.read_csv(filename, sep = '\t', usecols=['Query_name', 'Barcode']) for filename in inpathlist), ignore_index = True)

assert len(df) == len(df['Query_name'].unique()) # Read names are unique after concatenating dataframes

# Make a dataframe with the barcodes and a list of read ids corresponding to those barcodes
df_list = (df
           .groupby('Barcode')['Query_name'].apply(list).reset_index() # Group reads per barcode
           .assign(Counts = lambda df_: df_['Query_name'].str.len()) # Count number of reads per unique bc
          )


# Plot the histogram
# TODO: Clean it up!
max_counts = df_list['Counts'].max()
vmin, vmax = (0, int(np.percentile(df_list['Counts'], 99.9)))
nbins = max_counts if max_counts < 10 else min(200, vmax - vmin)
df_list['Counts'].hist(range=(vmin, vmax), bins = nbins)
plt.yscale("log")
plt.axvline(snakemake.params["threshold"], color="salmon")
plt.savefig(snakemake.output["barcode_hist"], dpi=150)

# Threshold barcodes by depth
df_list = (df_list
           .loc[df_list['Counts'] >= snakemake.params["threshold"]]
           .reset_index(drop=True) # Reset index after filtering
           .reset_index() # Reset index to generate the Barcode_id
           .rename(columns={'index': 'Barcode_id'})
          )

# Save 
df_list.to_csv(snakemake.output["inv_barcode_codebook"], sep='\t', index=False)