import numpy as np
import pandas as pd
import re
import subprocess
import multiprocessing
import numpy as np

from Bio.Seq import Seq
from Bio import SeqIO

SEED_LENGTH = 9
SPACER_LENGTH = 20
BIN_DICT = {
    0:[-0.5, 0.1], # Lower bound 
    1:[0.1, 0.2],
    2:[0.2, 0.3],
    3:[0.3, 0.4],
    4:[0.4, 0.5],
    5:[0.5, 0.6],
    6:[0.6, 0.7],
    7:[0.7, 0.8],
    8:[0.8, 0.9],
    9:[0.9, 1.5] # Cannot set upper bound bc there are fitness values >= 1
}


k = 2 # Number of gRNAs per bin

def get_bin_ids(df_: pd.DataFrame,
                col: str='relative fitness (mean)',
                bin_dict: dict=BIN_DICT) -> np.ndarray:
    '''
    Returns a list of bin indices for each row in the dataframe.
    The bins are equally spaced according to bin_dict
    '''

    fitness = df_[col]

    # Use list comprehension to get the conditions
    # V2
    conditions = [(fitness >= bin_dict[i][0]) & (fitness < bin_dict[i][1]) for i in range(len(bin_dict))]
    
    # V1 
    # conditions = [
    #     (fitness >= bin_dict[i][0]) & (fitness < bin_dict[i][1])
    #     for i in range(len(bin_dict) - 1)
    # ] + [
    #     fitness >= bin_dict[len(bin_dict) - 1]
    # ]

    # bin_indices = [0, 1, 2, 3, 4]
    bin_indices = list(bin_dict.keys())

    bins = np.select(conditions, bin_indices, default=-1)
    return bins

def generate_sample_single_gene(df_,
                                bin_dict=BIN_DICT,
                                k=2,
                                col='relative fitness (mean)',
                                bin_index_name='bin_index_fitness'):
    '''
    Subsamples the dataframe such that there are k gRNAs per bin for each gene.
    If there are less than k gRNAs in a given bin, it will select the closest ones to the middle of the bin
    until reaching k gRNAs.
    '''
    df_counts = df_.groupby([bin_index_name]).size()

    for i in list(bin_dict.keys()): # Suubstitute with bin index
        if i not in df_counts.index: # If there are no gRNAs in a bin, add it with count 0
            df_counts[i] = 0 

    df_counts = df_counts.sort_index()

    # Select random sample from subset of df_ such that df_counts is greater than or equal to k
    df_sample = (df_.loc[df_[bin_index_name]
                    .isin(df_counts[df_counts >= k].index)]
                    .groupby([bin_index_name])
                    .sample(n=k, replace=False)
                )

    # Get the remaining ones i.e. the ones whose index is not in df_sample
    df_remaining = (df_
                    .loc[~df_.index.isin(df_sample.index)]
                    .assign(chosen = lambda df_: df_[bin_index_name].isin(df_counts[df_counts < k].index)) # Mark the ones in the bins as chosen
    )

    # Only do this if there are remaining guides to select
    count_remaining = k - df_counts[df_counts < k]
    if count_remaining.sum() > 0:
        for i in count_remaining.index:
            n = count_remaining[i]
            middle = (bin_dict[i][0] + bin_dict[i][1]) / 2
            # print('i ={0}, n = {1}, middle = {2}'.format(i, n, middle))
            df_remaining = df_remaining.assign(distances = lambda df_: np.abs(df_[col] - middle))
            indices_closest = df_remaining[~df_remaining['chosen']].nsmallest(n, 'distances').index
            df_remaining.loc[indices_closest, 'chosen'] = True
            # print(df_remaining.loc[df_remaining['chosen'], ['gene', 'relative fitness (mean)']])
        df_remaining = df_remaining.drop(columns=['distances'])

        # Set the ones with the smallest distance that has not been chosen to chosen
        df_sample = pd.concat([df_sample, df_remaining.loc[df_remaining['chosen']]]).drop(columns=['chosen'])
        
    return df_sample

def check_restriction_site(in_seq: str,
                           rest_seq: str = "GGTCTC") -> bool:
    '''
    Check if a sequence (or its reverse complement) contains a restriction site
    Parameters:
        in_seq: str - input sequence
        rest_seq: str - restriction site sequence (default: BsaI site)
    '''
    test_str = Seq(in_seq.upper())
    rc_test_str = test_str.reverse_complement()
    if (rest_seq in test_str) or (rest_seq in rc_test_str):
        return True
    else:
        return False
    
def check_all_sites(in_seq, rest_seqs=["GAAGAC", "GGTCTC", "CGATCG"]):
    '''
    Check if a sequence contains any of the restriction sites
    Parameters:
        in_seq: str - input sequence
        rest_seqs: list - list of restriction sites (default: BbsI, BsaI, PvuI)
    
    '''
    out_test = False
    for rest_seq in rest_seqs:
        rest_site = check_restriction_site(in_seq, rest_seq=rest_seq)
        out_test = out_test or rest_site
    return out_test

def add_adapter_seqs(spacer: str) -> str:
    '''
    Add adapter sequences to the spacer for cloning
    Parameters:
        spacer: str - spacer sequence
    '''
    adapter_5p = 'CGTACGACGGAAGACATTCTCCTGTAATGTAAAAAATGGGCTCGTGTTGTACAATAAATGT'
    adapter_3p = 'GTTTTCGTCTTCTTAAGGTGCC'
    output_seq = adapter_5p + spacer + adapter_3p
    return output_seq

def find_seed_start_positions(query_seed, genome_seq):
    """
    Find the start positions of the seed sequence in the genome sequence.

    Parameters
    ----------
    query_seed : str
        The seed sequence to search for in the genome sequence.
    genome_seq : str
        The genome sequence to search in.
    """
    seed_start = []
    for m in re.finditer(query_seed+r'.GG', genome_seq):
        seed_start.append(m.start() + 1) # +1 for 0-based index
    
    return seed_start

def find_all_hits_seeds(df_seeds, upper_seq, lower_seq):
        """
        Find all the hits of the seeds in the genome sequence

        Parameters
        ----------
        df_seeds : pd.DataFrame
            DataFrame with the seeds to find in the genome sequence
        upper_seq : str
                Upper strand of the genome sequence
        lower_seq : str
                Lower strand of the genome sequence
        """
        assert(len(upper_seq) == len(lower_seq))
        GENOME_SIZE = len(upper_seq)

        df_seeds_up = (df_seeds
        .assign(seed_start = lambda df_: df_['seed'].apply(find_seed_start_positions, args=(upper_seq,)))
                #  start_pos_lower = lambda df_: df_['seed'].apply(find_seed_start_positions, args=(lower_seq,)))
        # Expand the data frame into a row for each seed start position, and a new column for the strand
        .explode('seed_start')
        .loc[lambda df_:~df_['seed_start'].isna()]
        .assign(matches_up_strand = True)   
        )

        df_seeds_lo = (df_seeds
                .assign(seed_start = lambda df_: df_['seed'].apply(find_seed_start_positions, args=(lower_seq,)))
                        #  start_pos_lower = lambda df_: df_['seed'].apply(find_seed_start_positions, args=(lower_seq,)))
                # Expand the data frame into a row for each seed start position, and a new column for the strand
                .explode('seed_start')
                .loc[lambda df_:~df_['seed_start'].isna()]
                .assign(matches_up_strand = False)
                .assign(seed_start = lambda df_: GENOME_SIZE - df_['seed_start'] + 1) # Desi
                )

        df_seeds_all = (pd.concat([df_seeds_up, df_seeds_lo], axis=0)
                        .sort_index()
                        .assign(seed_end = lambda df_: np.where(df_["matches_up_strand"],
                                                        df_["seed_start"] + SEED_LENGTH - 1, 
                                                        df_["seed_start"] - (SEED_LENGTH - 1)))
                        .assign(spacer_start = lambda df_: np.where(df_["matches_up_strand"],
                                                        df_["seed_start"] - (SPACER_LENGTH - SEED_LENGTH),
                                                        df_["seed_start"] + (SPACER_LENGTH - SEED_LENGTH)))
                        .assign(spacer_seq = lambda df_: df_.apply(lambda row: upper_seq[row["spacer_start"] - 1 : 
                                                                        row["seed_end"]]
                                                                if row["matches_up_strand"] 
                                                                else lower_seq[GENOME_SIZE - row['spacer_start'] : 
                                                                            GENOME_SIZE - (row['seed_end'] - 1)], axis=1))
                        # Prepare for bioframe
                        .assign(start = lambda df_: np.where(df_["matches_up_strand"],
                                                            df_["spacer_start"], 
                                                            df_["seed_end"]),
                                end = lambda df_  : np.where(df_["matches_up_strand"],
                                                            df_["seed_end"], 
                                                            df_["spacer_start"]),
                                chrom = '0')
                        .astype({"start": 'int64', "end": 'int64'})
        )

        del df_seeds_lo, df_seeds_up
        return df_seeds_all

import multiprocessing

# Define a function to be executed in parallel
def process_chunk(chunk_and_params):
    """
        Find all the hits of the seeds in the genome sequence

        Parameters
        ----------
        chunk_and_params : tuple
            Tuple with the chunk of the dataframe and the upper and lower strands of the genome sequence
    """
    chunk = chunk_and_params[0]
    upper_seq, lower_seq = chunk_and_params[1], chunk_and_params[2]
    
    return find_all_hits_seeds(chunk, upper_seq, lower_seq)

def find_all_hits_seeds_parallel(df_seeds, upper_seq, lower_seq):
    """
        Find all the hits of the seeds in the genome sequence

        Parameters
        ----------
        df_seeds : pd.DataFrame
            DataFrame with the seeds to find in the genome sequence
        upper_seq : str
                Upper strand of the genome sequence
        lower_seq : str
                Lower strand of the genome seque
    """

    # Define the number of processes to use
    num_processes = multiprocessing.cpu_count()

    # Split the dataframe into chunks for parallel processing
    chunks = np.array_split(df_seeds, num_processes)

    # Combine chunks and parameters as a list of tuples
    chunks_and_params = [(chunk, upper_seq, lower_seq) for chunk in chunks]

    # Create a pool of processes
    pool = multiprocessing.Pool(processes=num_processes)

    # Apply the function to each chunk in parallel
    results = pool.map(process_chunk, chunks_and_params)

    # Concatenate the results
    df_seeds_all = pd.concat(results)

    # Close the pool of processes
    pool.close()
    pool.join()

    return df_seeds_all.sort_values(by=['seed', 'seed_start'])

def subsample_spacers_per_bin(df: pd.DataFrame,
                              k: int = 2,
                              sort_by: list = ['n_off_targets', 'index_mismatch', 'binds_non_template_strand', 'n_nts_mismatched'],
                              sort_ascending: list = [True, True, True, False],
                              bin_index_name: str = 'bin_index_fitness'):
        '''
        Subsample spacers per bin
        1. Sort the dataframe by locus_tag, bin_index_name, and sort_by
        2. Select the top k spacers per bin that haven't been selected
        3. Assign the bin_index_name to the selected spacers
        4. Return the dataframe

        Parameters
        ----------
        df : pd.DataFrame
        Dataframe with the spacers
        k : int
        Number of spacers to select per bin
        sort_by : list
        List of columns to sort by
        sort_ascending : list
        List of booleans to indicate if the sorting is ascending or descending
        bin_index_name : str
        Name of the column with the bin index
        '''
        df = df.sort_values(
                by = ['locus_tag', bin_index_name] + sort_by,
                ascending = [True, True] + sort_ascending)
        # Among the ones that haven't been selected, select the k best ones per bin
        indices_chosen = df.loc[lambda df_: ~df_['selected']].groupby(['locus_tag', bin_index_name]).head(k).index
        return (df
                # indices_chosen = df.loc[lambda df_: ~df_['selected']].head(n_left_to_sample).index # HEREEE
                # .assign(selected = lambda df_: df_.index.isin(df_.groupby(['locus_tag', bin_index_name]).head(k).index)) # Select the top k from the existing bins
                .assign(selected = lambda df_: np.where(df_.index.isin(indices_chosen), True, df_['selected']))
                # Among the ones that are selected, assign the value of a column called bin_assigned to be equal to the bin_index_name
                .assign(bin_assigned = lambda df_: np.where(df_.index.isin(indices_chosen), df_[bin_index_name], None))
                .rename(columns={'bin_assigned': bin_index_name+'_assigned'})
        ) 

def generate_bin_sample_single_gene(df,
                                    bin_dict = BIN_DICT,
                                    k = 1,
                                    sort_by = ['n_off_targets', 'binds_non_template_strand_off_target', 'n_nts_mismatched',
                                                'agg_pam', 'n_t_in_start', 'offset_gt_1000'],
                                    sort_ascending = [True, True, False,
                                                        True, True, True],
                                    bin_index_name = 'bin_index_fitness'):
    '''
    Generate a sample of spacers per bin for a single gene using the closest distance to the bin index

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe with the spacers
    bin_dict : dict
        Dictionary with the bins
    k : int
        Number of spacers to select per bin
    sort_by : list
        List of columns to sort by
    sort_ascending : list
        List of booleans to indicate if the sorting is ascending or descending
    bin_index_name : str
        Name of the column with the bin index
    '''

    df_counts_per_bin = df.groupby([bin_index_name]).size()

    for i in list(bin_dict.keys()): # Suubstitute with bin index
        if i not in df_counts_per_bin.index: # If there are no gRNAs in a bin, add it with count 0
            df_counts_per_bin[i] = 0 

    df_counts_per_bin = df_counts_per_bin

    df_counts_per_bin = (
        df_counts_per_bin
            .sort_index()
            .to_frame()
            .reset_index()
            .rename(columns={0: 'n_total'})
            .assign(n_left_to_sample = lambda df_: np.where(df_['n_total'] >= k, 0, k - df_['n_total']))
    )

    for _, row in df_counts_per_bin.iterrows():
        bin_index = row[bin_index_name]
        n_left_to_sample = row['n_left_to_sample']
        if n_left_to_sample == 0:
            continue
        # print(bin_index,n_left_to_sample)
        # Get the closest bin that still has gRNAs
        df = (df
                .assign(distance = lambda df_: abs(df_[bin_index_name] - bin_index))
                .sort_values(by = ['distance'] + sort_by,
                             ascending = [True] + sort_ascending)
                            # by = ['distance', 'n_off_targets', 'binds_non_template_strand', 'n_nts_mismatched'], 
                            # ascending = [True, True, True, False])
        )

        # Set the n_left_to_sample first ones for the minimum distance to chosen
        indices_chosen = df.loc[lambda df_: ~df_['selected']].head(n_left_to_sample).index
        df = (df
            .assign(selected = lambda df_: np.where(df_.index.isin(indices_chosen), True, df_['selected']))
            .assign(**{bin_index_name+'_assigned': lambda df_: np.where(df_.index.isin(indices_chosen), bin_index, df_[bin_index_name+'_assigned'])})
        )
        # df.loc[indices_chosen, 'selected'] = True
        # df.loc[indices_chosen
        # , 'bin_assigned'] = bin_index
        # break
    # If there is a column named distance, drop it
    if 'distance' in df.columns:
        df = df.drop(columns = ['distance'])
    df.sort_values(by = [bin_index_name] + sort_by,
                     ascending = [True] + sort_ascending)
        # by = [bin_index_name, 'n_off_targets', 'binds_non_template_strand', 'n_nts_mismatched'],
                #    ascending = [True, True, True, False])
    return df

def count_differences(seq1, seq2):
    """
    Count the number of differences between two sequences of the same length.
    Returns nan if either one of them is nan
    """
    if pd.isna(seq1) or pd.isna(seq2):
        return np.nan
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def run_grna_finding_pipeline(locus_tag_list, 
                              path_to_script,
                              path_to_genbank,
                              path_to_params,
                              path_log,
                              path_to_output,
                              get_mismatches = False):
    """
    Run the gRNA finding pipeline using code from the paper
    locus_tag_list: list of locus tags
    
    path_to_script: path to the script
    path_to_genbank: path to the genbank file
    path_to_params: path to the model parameters
    path_log: path to the log file
    path_to_output: path to the output file
    get_mismatches: whether to get the mismatches or not
    """


    # Find
    with open(path_log, 'w') as f:
        subprocess.run(['python', path_to_script, 
                        '--genbank', path_to_genbank,
                        '--locus_tag', ','.join(locus_tag_list),
                        '--model_param', path_to_params, 
                        '--step', 'find',
                        '--file_find', path_to_output+'_find.csv'
                        ], stdout=f, stderr=f)
    # Filter
    with open(path_log, 'w') as f:
        subprocess.run(['python', path_to_script, 
                        '--genbank', path_to_genbank,
                        '--locus_tag', ','.join(locus_tag_list),
                        '--model_param', path_to_params, 
                        '--step', 'filter',
                        # '--sgrna_remove', ,
                        # '--pam_remove', 'ANN',
                        # '--gc_upper', '0.6',
                        '--off_target_upper', '10',
                        '--offset_upper', '1000000', # Infinite
                        '--file_find', path_to_output+'_find.csv',
                        '--file_filter', path_to_output+'_filter.csv'
                        ], stdout=f, stderr=f)
    if get_mismatches:
        # Mismatch
        with open(path_log, 'w') as f:
            subprocess.run(['python', path_to_script, 
                            '--genbank', path_to_genbank,
                            '--locus_tag', ','.join(locus_tag_list),
                            '--model_param', path_to_params, 
                            '--step', 'mismatch',
                            '--file_find', path_to_output+'_find.csv',
                            '--file_filter', path_to_output+'_filter.csv',
                            '--file_mismatch', path_to_output+'_mismatch.csv'
                            # '--sgrna_remove', ,
                            ], stdout=f, stderr=f)
            

def load_and_process_dfs_from_pipeline(path_to_output,
                                       df_genes,
                                       bin_dict=BIN_DICT):
    """
    Load the filtered (fully on-target) and the mismatch dataframes that the pipeline generates.
    Process the dataframes (read recipe) and returns a dataframe with all the spacers and information about them

    Parameters:
    path_to_output (str): Path to the output files of the pipeline
    df_genes (pd.DataFrame): DataFrame with the information of all genes
    bin_dict (dict): Dictionary with the bins
    """

    # Load the filtered (fully on-target) and the mismatch dataframes
    df_f = (pd.read_csv(path_to_output+'_filter.csv',
                    usecols=['locus_tag', 'gene_name', 'strand', 'n', 'offset',
                            'position_offset', 'gc', 'sgrna', 'pam',
                            'downstream', 'sequence', 'seed', 'off_target_n'])
        .rename(columns={'sgrna': 'full_sgrna', 'gene_name': 'gene'})
        .astype({'strand': 'category',
                'n': 'int32',
                'offset': 'int32',
                'position_offset': 'int32',
                'pam': 'category',
                'off_target_n': 'int32'})
        .assign(predicted_efficacy = 1,
                index_mismatch = 0,
                mismatch_sgrna = lambda df_: df_['full_sgrna'])
        
        .assign(bin_index = lambda df_:get_bin_ids(df_,
                                                    col='predicted_efficacy',
                                                    bin_dict=bin_dict)) # Assign bin indices to each gRNA
        .astype({'predicted_efficacy': 'float64', 'bin_index': 'int8'})                                                                   
    )

    df_m = (pd.read_csv(path_to_output+'_mismatch.csv',
                    usecols=['locus_tag', 'n', 'position', 'original', 'substitution', 'full_sgrna', 'mismatch_sgrna', 'predicted efficacy'])
        .rename(columns={'predicted efficacy': 'predicted_efficacy'})
        .astype({'original': 'category', 'substitution': 'category'})
        
        .assign(index_mismatch = lambda df_: df_['n'].str.extract(r'\((\d+)\)').astype('int16'))
        .assign(n = lambda df_: df_['n'].str.extract(r'(\d+)\(').astype('int32'))
    #       .merge(df_aux[['locus_tag', 'gene']], on='locus_tag', how='left') # Recover gene names
        # .merge(df_genes.loc[df_genes['category_gene']=='auxotroph',['locus_tag', 'gene']], on='locus_tag', how='left')
        .assign(bin_index = lambda df_:get_bin_ids(df_,
                                                    col='predicted_efficacy',
                                                    bin_dict=bin_dict)) # Assign bin indices to each gRNA
        .astype({'bin_index': 'int8'})

        )

    # Concatenate and merge both dfs
    df = (pd.concat([df_m, df_f], axis=0)
    # Drop the columns that wil be superseded by the value of the filtered df (no mismatches)
    .drop(columns = ['strand', 'offset', 'position_offset', 'gc', 'pam', 'downstream', 'sequence', 'seed', 'off_target_n'])
    .merge((df_f
            # Drop the columns that will be superseded by the value of the mismatch df
            .drop(columns = ['gene', 'full_sgrna', 'predicted_efficacy', 'index_mismatch',
                            'bin_index', 'mismatch_sgrna']
            )),
            on=['locus_tag', 'n'])

    .loc[lambda df_: df_['mismatch_sgrna'].duplicated(keep=False) == False] # Remove duplicate sgRNAs (both of them)    
    .loc[lambda df_: df_['mismatch_sgrna'].apply(check_all_sites) == False] # Remove rows with restriction sites
    .loc[lambda df_: df_['mismatch_sgrna'].apply(add_adapter_seqs).apply(check_all_sites, args=[["GGTCTC", "CGATCG"]]) == False] # Remove rows with BsaI and PvuI sites
    .loc[lambda df_: df_.bin_index != -1] # Remove gRNAs that were not assigned to a 
    .sort_values(by=['locus_tag', 'n', 'index_mismatch'])
    .reset_index(drop=True)

    .rename(columns={'mismatch_sgrna': 'query_spacer'})
    .assign(query_end_position = lambda df_: df_['position_offset'])
    
    # CORRECT LOCUS_TAGS AND GENE NAMES
    .rename(columns={'locus_tag': 'locus_tag_file'})
    .drop(columns = ['gene'])
    .merge(df_genes.loc[:,['locus_tag_file', 'locus_tag', 'gene']],
            left_on='locus_tag_file',
            right_on='locus_tag_file',
            how='left')
    .drop(columns= ['locus_tag_file'])
    .loc[:, lambda df_: ['locus_tag', 'gene']+ [col for col in df_.columns if col not in ['locus_tag', 'gene']]]

    # Apply booleans for filtering
    # Count the number of gs and cs on the 7nth and 8th position of the downstream sequence
    .assign(gc_downstream_8_9 = lambda df_: df_['downstream'].str[7:9].str.count('G') + df_['downstream'].str[7:9].str.count('C'),
            offset_gt_1000 = lambda df_: df_['offset'] > 1000,
            # Assign a column that counts the number of Ts in the first 3 bases of the query_spacer
            n_t_in_start = lambda df_: df_['query_spacer'].str[:3].str.count('T'),
            agg_pam = lambda df_: df_['pam'].str.contains('AGG'))

    )
    return df

def change_commas_to_semicolon_in_locus_tag(input_gb, output_gb):
    """
    Change commas to semicolons in the locus_tag and label qualifiers of a GenBank file.

    Args:
        input_gb (str): The input GenBank file.
        output_gb (str): The output GenBank file.
    """
    # Open the GenBank file
    with open(input_gb, 'r') as input_file:
        records = SeqIO.parse(input_file, 'genbank')

        # Create a new list to store the modified records
        modified_records = []

        # Iterate over each record in the GenBank file
        for record in records:
            # Iterate over each feature in the record
            for feature in record.features:
                # Check if the feature has a locus tag and replace commas
                if 'locus_tag' in feature.qualifiers:
                    feature.qualifiers['locus_tag'] = [tag.replace(',', ';') for tag in feature.qualifiers['locus_tag']]

                # Check if the feature has a label and replace commas
                if 'label' in feature.qualifiers:
                    feature.qualifiers['label'] = [label.replace(',', ';') for label in feature.qualifiers['label']]

            # Add the modified record to the list
            modified_records.append(record)

    # Write the modified records to a new GenBank file
    with open(output_gb, 'w') as output_file:
        SeqIO.write(modified_records, output_file, 'genbank')

def load_genes(path_genbank):
    """
    Load genes from a genbank file into a dataframe

    Args:
    path_genbank (str): Path to the genbank file
    """
    rec= SeqIO.read(path_genbank, "genbank")

    upper_seq = str(rec.seq)
    lower_seq = str(rec.seq.reverse_complement())
    GENOME_SIZE = len(upper_seq)

    # Load features in a dictionary
    if not rec.features:
        raise ValueError("No features found in the genbank file")

    dict_genes = {"locus_tag_file":[],
                "gene":[],
                "strand":[],
                "start":[],
                "end":[],
                #   "sequence":[]
                }

    for feature in rec.features:
        
        if feature.type == "gene" and "locus_tag" in feature.qualifiers:
            # print(i)
            # i+=1
            # dict_genes["locus_tag_file"].append(feature.qualifiers.get("locus_tag")[0])
            dict_genes["locus_tag_file"].append(feature.qualifiers.get("locus_tag", ["NA"])[0])
            
            if feature.qualifiers.get("label"):
                dict_genes["gene"].append(feature.qualifiers.get("label")[0])
            else:
                dict_genes["gene"].append(feature.qualifiers.get("locus_tag", ["NA"])[0])
                # dict_genes["locus_tag_file"].append(feature.qualifiers.get("locus_tag", ["NA"])[0])

            # dict_genes["gene"].append(feature.qualifiers.get("label")[0])
            dict_genes["strand"].append(feature.location.strand)
            dict_genes["start"].append(int(feature.location.start))
            dict_genes["end"].append(int(feature.location.end))
            # dict_genes["sequence"].append(str(feature.extract(rec.seq)))
    
    # Create dataframe and fiter
    df_genes = (pd.DataFrame(dict_genes)
                .assign(length = lambda df_: df_["end"] - df_["start"],
                        locus_tag = lambda df_: df_["locus_tag_file"].str.split(';').str[0])
            #     .assign(locus_tag_right = lambda df_: df_["locus_tag_file"].str.split(',').str[1])
                .loc[lambda df_:~df_['locus_tag'].str.contains('new')]
                .loc[lambda df_:~df_['locus_tag'].str.contains('misc')]
                .assign(locus_tag = lambda df_: df_['locus_tag'].apply(lambda x: re.sub(r'U(\d)', r'U_\1', x)))
                # Relocate columns
                .loc[:, lambda df_:['locus_tag', 'gene'] + [col for col in df_.columns if col not in ['locus_tag','gene']]]
                .assign(matches_up_strand = lambda df_:df_['strand'] == 1)
                .astype({'matches_up_strand': 'boolean'})
    )
                # .to_csv(headpath+'genes_test.csv', index=False)
    return df_genes


def classify_genes(df_genes,
                   headpath,
                   filename_essentials,
                   filename_auxotrophs):
    """
    Classify genes into essential, auxotroph, r_rna, t_rna, and nonessential

    Args:
    df_genes: DataFrame
    headpath: (str) Path to the directory where the files are located
    filename_essentials: (str) Name of the file containing essential genes (from paper)
    filename_auxotrophs: (str) Name of the file containing auxotroph genes (from paper)
    """
    genes_ess = (pd.read_excel(headpath + filename_essentials,
                        sheet_name="relative fitness (bsu data)",
                        usecols=['locus_tag', 'gene'])
                    .rename(columns={'gene':'gene_name_data'})
                    .drop_duplicates(subset='locus_tag', keep='first')
                    .assign(locus_tag = lambda df_: df_['locus_tag'].str.slice(0, 3)+'_'+df_['locus_tag'].str.slice(3)) # Add underscore to locus tag
    )

    # Read Excel file
    # table_url = "https://ars.els-cdn.com/content/image/1-s2.0-S2405471220303665-mmc4.xlsx"
    genes_aux = (pd.read_excel(headpath+filename_auxotrophs,
                        sheet_name="D",
                        skiprows=16,
                        usecols = ['Locus tag', 'gene']
                        )
            # Add _ to the Locus tag column
            .rename(columns={"Locus tag": "locus_tag",
                            "gene": "gene_name_data"})
            .assign(locus_tag = lambda df_: df_['locus_tag'].str.slice(0, 3)+'_'+df_['locus_tag'].str.slice(3))
            .loc[lambda df_:~df_['locus_tag'].isin(genes_ess['locus_tag'])]
        ) 

    df_genes = (df_genes
                .merge(genes_ess, on='locus_tag', how='left')
                .assign(category_gene = lambda df_:np.where(df_['locus_tag'].isin(genes_ess['locus_tag']),
                                                            'essential_from_dataset',
                                                            pd.NA))
                .merge(genes_aux, on='locus_tag', how='left')
                .assign(category_gene = lambda df_:np.where(df_['locus_tag'].isin(genes_aux['locus_tag']),
                                                            'auxotroph',
                                                            df_['category_gene']))
                .assign(gene_name_data = lambda df_: df_['gene_name_data_x'].combine_first(df_['gene_name_data_y']))
                .drop(columns = ['gene_name_data_x', 'gene_name_data_y'])
                # .loc[lambda df_:df_['category_gene'].isin(['auxotroph', 'essential_from_dataset'])]

                .assign(category_gene = lambda df_: np.where(df_['locus_tag'].str.contains(r'BSU_[r]'),
                                                'r_rna',
                                                df_['category_gene']))
                .assign(category_gene = lambda df_: np.where(df_['locus_tag'].str.contains(r'BSU_[t]'),
                                                            't_rna',
                                                            df_['category_gene']))
                                                            #   .loc[lambda df_:df_['category_gene'] == 't_rna']
                .assign(category_gene = lambda df_: np.where(df_['category_gene'].isna(),
                                                            'nonessential',
                                                            df_['category_gene']))
    )

    return df_genes