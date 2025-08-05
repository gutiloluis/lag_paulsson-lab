import pandas as pd
import glob
import re
import os
import numpy as np
import skimage as sk


def parse_and_save_futures_results(futures_dict, function_name, headpath):
    """
    Parse the results of futures and save them to a DataFrame and a pickle file.

    Args:
        futures_dict (dict): Dictionary with indices as keys and futures as values.
        function_name (str): Name of the function used to generate the futures.
        headpath (str): Base path to save the pickle file.

    Returns:
        pd.DataFrame: DataFrame containing the status, done, result, and exception of each future.
    """
    def parse_futures_dict_entry_to_row(future):
        """
        Parse a single future entry to a dictionary row.

        Args:
            future (dask.distributed.Future): Future object to parse.

        Returns:
            dict: Dictionary containing the status, done, result, and exception of the future.
        """
        try:
            result = future.result()
        except Exception as e:
            result = str(e)
        future_dict_row = {
            'status': future.status,
            'done': future.done(),
            'result': result,
            'exception': str(future.exception()),
            'key': future.key
        }
        return future_dict_row

    futures_dict_to_table = {
        index: parse_futures_dict_entry_to_row(future)
        for index, future in futures_dict.items()
    }

    df_futures = (
        pd.DataFrame
        .from_dict(futures_dict_to_table, orient='index')
        .astype({
            'status': 'category',
            'done': bool,
            'result': 'object',
            'exception': 'object'
        })
        .reset_index()
        .rename(columns={'index': 'index'})
        .sort_values(by='index', ascending=True)
    )

    df_futures.to_pickle(f"{headpath}/futures_{function_name}.pkl")

    return df_futures


def deploy_function_in_dask(dask_controller, headpath, function, indices, *args, **kwargs):
    """
    Deploy a function in Dask for each index in the list.

    Args:
        dask_controller: Object containing a daskclient attribute (e.g., from distributed.Client).
        headpath (str): Base path to save the result pickle file.
        function (callable): Function to deploy.
        indices (list): List of indices to process.
        *args: Additional positional arguments for the function.
        **kwargs: Additional keyword arguments for the function.

    Returns:
        pd.DataFrame: Results DataFrame with future metadata and results.
    """
    futures_dict = {}
    for index in indices:
        future = dask_controller.daskclient.submit(
            function,
            index,
            *args,
            **kwargs
        )
        futures_dict[index] = future

    df_futures = parse_and_save_futures_results(futures_dict, function.__name__, headpath)
    return df_futures

def group_filenames_per_fov(headpath, pattern="*BF_*"):
    """
    Groups filenames in the specified directory by their FOV number.
    Args:
        headpath (str): The directory path where the files are located.
        pattern (str): The pattern to match filenames. Default is "*BF_*".
    Returns:
        dict: A sorted dictionary where keys are FOV numbers and values are lists of filenames.
    """
    filenames = glob.glob(headpath + pattern)
    filenames_per_fov = {}
    fov_pattern = re.compile(r"fov=(\d+)")

    for filename in filenames:
        match = fov_pattern.search(os.path.basename(filename))
        if match:
            fov_number = int(match.group(1))
            if fov_number not in filenames_per_fov:
                filenames_per_fov[fov_number] = []
            filenames_per_fov[fov_number].append(filename)
    # Also sort keys
    filenames_per_fov = dict(sorted(filenames_per_fov.items()))
    sorted_filenames_per_fov = {}
    # Sort filenames
    for fov_number, filenames in filenames_per_fov.items():
        sorted_filenames = sorted(filenames, key=lambda x: int(re.search(r'_t=(\d+)', x).group(1)))
        sorted_filenames_per_fov[fov_number] = sorted_filenames
    return sorted_filenames_per_fov

def generate_output_filenames_dict_per_fov(filenames_per_fov_in, headpath_out):
    filenames_per_fov_out = {fov: [] for fov in filenames_per_fov_in.keys()}
    for fov, filenames in filenames_per_fov_in.items():
        for filename in filenames:
            filename_out = os.path.join(headpath_out, os.path.basename(filename))
            filenames_per_fov_out[fov].append(filename_out)

    return filenames_per_fov_out

def group_filenames_per_fov_channel_all_timepoints(headpath, pattern="*fov*", n_timepoints= 10, channels = ['BF', 'RFP', 'Cy5', 'Cy7']):
    '''
    Groups filenames per FOV and channel for all timepoints.
    Args:
        headpath (str): Path to the directory containing the files.
        pattern (str): Pattern to match the filenames.
        n_timepoints (int): Number of timepoints expected per channel.
    Returns:
        dict: A dictionary where keys are FOVs and values are dictionaries with channels as keys and lists of filenames as values.
    '''
    filenames_per_fov_all_channels_all_timepoints = group_filenames_per_fov(headpath=headpath, pattern=pattern)
    filenames_per_fov_channel_all_timepoints = {}
    for fov, filenames in filenames_per_fov_all_channels_all_timepoints.items():
        filenames_per_fov_channel_all_timepoints[fov] = {channel: [] for channel in channels}
        for filename in filenames:
            channel = re.search(r'config=([A-Za-z0-9]+)_', filename).group(1)
            assert channel in channels, f"Channel {channel} not in {channels}"
            filenames_per_fov_channel_all_timepoints[fov][channel].append(filename)
    assert len(filenames_per_fov_channel_all_timepoints) == len(filenames_per_fov_all_channels_all_timepoints)
    assert np.all(
        (np.array([len(filenames_per_fov_channel_all_timepoints[fov][channel]) for fov in filenames_per_fov_channel_all_timepoints for channel in channels]) == n_timepoints)
    )
    return filenames_per_fov_channel_all_timepoints

## These two functions are for the metadata files before extraction
def add_fake_channel_to_metadata(df_metadata, config_template='BF', config_new='BFSmoothed'):
    """
    Adds a new configuration to the metadata DataFrame for BF smoothed data.
    Args:
        df_metadata (pd.DataFrame): The original metadata DataFrame.
        config_template (str): The configuration to use as a template for the new config.
        config_new (str): The name of the new configuration to be added.
    Returns:
        pd.DataFrame: The updated metadata DataFrame with the new configuration added.
    """
    df_subset_new = (\
        df_metadata
        .loc[df_metadata['config'] == config_template, :]
        .assign(config=config_new)
    )
    return (pd
            .concat([df_metadata, df_subset_new], ignore_index=True)
            .sort_values(by=['fov', 't', 'config'], ignore_index=True)
    )

def add_fake_channel_to_all_metadata_files(metadata_filenames, config_template='BF', config_new='BFSmoothed'):
    """
    Adds a new configuration to all metadata files in the specified list.
    Args:
        metadata_filenames (list): List of metadata file paths.
    Returns:
        None
    """
    for metadata_filename in metadata_filenames:
        df_metadata = pd.read_hdf(metadata_filename)
        df_metadata = add_fake_channel_to_metadata(df_metadata, config_template, config_new)
        df_metadata.to_hdf(metadata_filename, key='data', mode='w')
    return 0


def correct_drift_single_image(image, cumulative_shift_coords_t):
    shift_y, shift_x = cumulative_shift_coords_t
    transform = sk.transform.AffineTransform(translation=(-shift_x, -shift_y))
    image_corrected = sk.transform.warp(image,
                                        transform,
                                        mode='constant',
                                        preserve_range=True,
                                        cval=np.median(image)).astype(image.dtype)
    return image_corrected