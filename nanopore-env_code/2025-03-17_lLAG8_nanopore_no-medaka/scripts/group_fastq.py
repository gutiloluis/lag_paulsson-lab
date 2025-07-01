import os
import shutil
import gzip

input_dir = snakemake.input["fastq_path"]
output_dir = snakemake.output["fastq_chunkpath"]
n = snakemake.params["fastq_chunksize"]

    
# List all .fastq.gz files in the input directory
fastq_files = [f for f in os.listdir(input_dir) if f.endswith('.fastq.gz')]

# Sort the files to maintain order
fastq_files.sort()

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Initialize chunk counter
chunk_index = 0

# Iterate over the files in chunks of n
for i in range(0, len(fastq_files), n):
    # Get the current chunk of files
    chunk_files = fastq_files[i:i + n]
    
    # Output file name
    output_filename = os.path.join(output_dir, f"chunk_{chunk_index}.fastq")
    
    # Open the output file for writing
    try:
        with open(output_filename, 'wb') as outfile:
            # Iterate over each file in the chunk
            for filename in chunk_files:
                filepath = os.path.join(input_dir, filename)
                # Open the .fastq.gz file for reading
                with gzip.open(filepath, 'rb') as infile:
                    # Read the content and write it to the output file
                    outfile.write(infile.read())
    except Exception as e:
        print(f"Error processing chunk {chunk_index}: {e}")
        continue
    
    # Increment the chunk counter
    chunk_index += 1

