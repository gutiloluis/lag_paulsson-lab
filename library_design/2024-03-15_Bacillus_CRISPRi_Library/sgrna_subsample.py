import numpy as np
import pandas as pd

# # bin dict must be in ascending order. Upper bound should be a number, not an interval
# bin_dict = {
#     0:[0.0, 0.2],
#     1:[0.2, 0.4],
#     2:[0.4, 0.6],
#     3:[0.6, 0.8],
#     4: 0.8 # Cannot set upper bound bc there are fitness values >= 1
# }

bin_dict = {
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

# bin_dict = {
#     0:[0.0, 0.1],
#     1:[0.1, 0.2],
#     2:[0.2, 0.3],
#     3:[0.3, 0.4],
#     4:[0.4, 0.5],
#     5:[0.5, 0.6],
#     6:[0.6, 0.7],
#     7:[0.7, 0.8],
#     8:[0.8, 0.9],
#     9: 0.9 # Cannot set upper bound bc there are fitness values >= 1
# }

k = 2 # Number of gRNAs per bin

def get_bin_ids(df_: pd.DataFrame,
                col: str='relative fitness (mean)',
                bin_dict: dict=bin_dict) -> np.ndarray:
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
                                bin_dict=bin_dict,
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