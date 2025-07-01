import os
import shutil

pod5path = snakemake.input["pod5path"]
pod5chunkpath = snakemake.output["pod5chunkpath"]
pod5_chunksize = snakemake.params["pod5_chunksize"]

if os.path.exists(pod5chunkpath):
    shutil.rmtree(pod5chunkpath)
os.makedirs(pod5chunkpath)

pod5_file_list = list(
    sorted([item for item in os.listdir(pod5path) if item[-5:] == ".pod5"])
)

for file_idx, item in enumerate(pod5_file_list):
    chunk_id = file_idx // pod5_chunksize
    pod5_chunkid_path = pod5chunkpath + "/chunk_" + str(chunk_id)
    if not os.path.exists(pod5_chunkid_path):
        os.makedirs(pod5_chunkid_path)
    shutil.copyfile(pod5path + "/" + item, pod5_chunkid_path + "/" + item)
