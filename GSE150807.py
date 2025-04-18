#&=================================================================================================*
#&  WRITTEN BY : Asmaa Yassin                                                                      *
#&  DATE       : April 17, 2025                                                                    *
#&-------------------------------------------------------------------------------------------------*
#& DESCRIPTION: This script loads and pre-processes raw RNA-seq data GSE150807                     *
#& removes missing values and genes with zero data and the cleansed data is exported to a CSV file.*
#& Metadata csv file is created containing the sample condition for later anylsis                  *
#&=================================================================================================*


import pandas as pd
import numpy as np


### Load Raw Count File ###
input_counts_file = "path/to/GSE150807_raw_counts_GRCh38.p13_NCBI.tsv"
df = pd.read_csv(input_counts_file, sep='\t')
print(df.head()) # Preview first few rows

df.set_index("GeneID", inplace=True) # Set GeneID column as the index
print(df.head()) # confirming new index
print(df.dtypes)  # Check if values are integer type

### Clean and Filter Data ###

print("Missing values per column:\n", df.isnull().sum()) # Check for missing values per column
df_clean = df.dropna() # Remove rows with any missing values (NaNs)
df_filtered = df_clean[(df_clean != 0).any(axis=1)] # Remove genes with all-zero counts across all samples

print("Original shape:", df.shape)  # check original data size
print("After dropping NAs:", df_clean.shape)  # check data size after removing (NaNs)
print("After removing all-zero rows):", df_filtered.shape) # check data size after removing all zeros

### Save Cleaned Data ###
output_counts_file = "path/to/GSE150807 raw counts filtered.csv" #set path to save output file
df_filtered.to_csv(output_counts_file)  # save cleaned data to csv
print(f"\n Filtered counts matrix saved to: {output_counts_file}")

### Create and Save Metadata File ###
metadata = pd.DataFrame({
    "sample": ["GSM4559030", "GSM4559031", "GSM4559032",
               "GSM4559033", "GSM4559034", "GSM4559035"],
    "condition": ["sensitive", "sensitive", "sensitive",
                  "resistant", "resistant", "resistant"],
    "replicate": [0, 1, 2, 0, 1, 2] })
metadata_output = "path/to/GSE150807_metadata.csv"
metadata.to_csv(metadata_output, index=False)
print(f"Metadata saved to: {metadata_output}")
