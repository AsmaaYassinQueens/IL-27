#&=================================================================================================*
#&  WRITTEN BY : Asmaa Yassin                                                                      *
#&  DATE       : April 17, 2025                                                                    *
#&-------------------------------------------------------------------------------------------------*
#& DESCRIPTION:This script merges multiple gene expression files from GSE140440 into one dataframe,*
#& renames columns, and checks mismatches with metadata to ensure all sample names match before    *
#& analysis.                                                                                       *
#&                                                                                                 *
#&=================================================================================================*


import os
import pandas as pd

directory = "path/to/GSE140440_RAW" # Set the directory containing the gene expression count files

### Listing and Filtering Count Files ###
files = os.listdir(directory)  # Retrieves the filenames in the directory
print(f"Found {len(files)} files in the directory.") #verifying # of files found
counts_files = [f for f in files if f.endswith(".counts") or f.endswith(".counts.txt")] # Filter to include only count files (with ".counts" or ".counts.txt")
print(f"Found {len(counts_files)} count files: {counts_files[:5]}")  #previewing first 5 files

### Reading and Merging Count Files ###
merged_df = None  #declaring a new variable to store the combined data
for file in counts_files: # Loop through each count file to read and merge
    sample_name = file.split(".")[0]  # Extract the sample name
    filepath = os.path.join(directory, file)
    temp_df = pd.read_csv(filepath, sep="\t", header=None, names=["Gene", sample_name])
    print(f"Reading file: {file}, Shape: {temp_df.shape}")  # Print file name and its dimensions
    if merged_df is None:
        merged_df = temp_df
    else:
        merged_df = merged_df.merge(temp_df, on="Gene", how="outer") # Merge on the "Gene" column using an outer join


### Saving Final Merged DataFrame to a CSV File ###
if merged_df is not None:
    output_file = os.path.join(directory, "GSE140440_merged_matrix.csv")  # Define output file path
    merged_df.to_csv(output_file, index=False)
    print(f"Merged file is saved as: {output_file}")  # Printing confirmation message
else:
    print("No data was merged, check file formats.")  #Printing a warning message if no data was merged

### set path to merged matrix and metadata ###
counts_file = "path/to/GSE140440_merged_matrix.csv"
metadata_file = "path/to/GSE140440_metadata.csv"

### Load and Clean Merged Data ###
df_counts = pd.read_csv(counts_file)
df_counts.columns = [col.split("_")[0] if "_" in col else col for col in df_counts.columns] # Rename columns by keeping only what's before the underscore
print("After renaming:", df_counts.columns[:5].tolist()) # Print first 5 column names after renaming
df_counts.to_csv(counts_file, index=False) # Save the file after renaming

df_metadata = pd.read_csv(metadata_file) # Load metadata

# Extract sample names from metadata and counts matrix
metadata_samples = set(df_metadata["Sample ID"])
count_samples = set(df_counts.columns[1:])  # Exclude 'Gene' column

# Finding mismatches between metadata and count matrix
missing_in_counts = metadata_samples - count_samples
missing_in_metadata = count_samples - metadata_samples

# Print mismatches
if missing_in_counts:
    print(f"Samples in metadata but missing from counts: {missing_in_counts}")
if missing_in_metadata:
    print(f"Samples in counts but missing from metadata: {missing_in_metadata}")
if not missing_in_counts and not missing_in_metadata:
    print("All sample names match between metadata and counts matrix and are ready for analysis!")
