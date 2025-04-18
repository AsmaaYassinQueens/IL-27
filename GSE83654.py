#&=================================================================================================*
#&  WRITTEN BY : Asmaa Yassin                                                                      *
#&  DATE       : April 17, 2025                                                                    *
#&-------------------------------------------------------------------------------------------------*
#& DESCRIPTION: This script merges multiple gene expression files from GSE83654,                   *
#& calculates correlation between samples, filters for significant genes,                          *
#& and identifies up and downregulated over time for PC3, DU145, and LNCaP cell lines              *
#& merged data. Correlation matrix, and upregulated and downregulated gene lists                   *
#& are returned as CSV files                                                                       *
#&=================================================================================================*

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def process_file(file_path):
    """
    Read and process single gene expression file:
      - Read columns 'GeneName' and 'LogRatio'
      - Drop rows with missing GeneName
      - Find the mean of LogRatio values for duplicate genes
    """
    data = pd.read_csv(
        file_path,
        delimiter='\t',
        skiprows=9,   #skip first 9 rows (metadata)
        usecols=['GeneName', 'LogRatio'],
        dtype={'GeneName': 'str', 'LogRatio': 'float32'}
    )
    data = data.dropna(subset=['GeneName']) # Drop missing gene name rows
    data = data.groupby('GeneName', as_index=False)['LogRatio'].mean() #find mean LogRatio for duplicate genes
    return data


def extract_label_from_filename(filename):
    """Extract a label from the file name by removing its extension."""
    base = os.path.basename(filename)
    label, _ = os.path.splitext(base)
    return label


def process_all_files_in_folder(folder_path):
    """
    1. Loop over all .txt files in `folder_path`.
    2. Process each file into a DataFrame.
    3. Rename 'LogRatio' column based on the filename.
    4. Merge all DataFrames on 'GeneName' (inner join).
    """
    txt_files = glob.glob(os.path.join(folder_path, '*.txt'))
    if not txt_files:
        raise FileNotFoundError(f"No .txt files found in {folder_path}")

    dataframes = []
    for file_path in txt_files:
        df = process_file(file_path)
        label = extract_label_from_filename(file_path)
        df.rename(columns={'LogRatio': f'_{label}'}, inplace=True)
        dataframes.append(df)

    # Merge all DataFrames on 'GeneName'
    merged_data = dataframes[0]
    for df in dataframes[1:]:
        merged_data = merged_data.merge(df, on='GeneName', how='inner', sort=False)

    return merged_data


if __name__ == "__main__":
    input_folder = "path/to/input_folder" # Set path of folder containing .txt files
    output_file = "path/to/GSE83654_data.csv"
    output_filtered = "path/to/GSE83654_filtered_genes.csv"
    output_up= "path/to/GSE83654_upregulated.csv"
    output_down = "path/to/GSE83654_downregulated.csv"

    ### Merge files ###
    merged_data = process_all_files_in_folder(input_folder)
    print("Merged Data:")
    print(merged_data.head())
    print("Number of genes in merged dataset:", len(merged_data))
    merged_data.to_csv(output_file, index=False)
    print(f"Merged data saved to {output_file}")


    ### Correlation Matrix ###
    logratio_cols = [col for col in merged_data.columns if col.startswith("_")] # Identify LogRatio columns (the ones starting with "_")

    if logratio_cols: # Display the correlation matrix
        corr_matrix = merged_data[logratio_cols].corr()
        print("\nCorrelation Matrix:")
        print(corr_matrix)
        corr_matrix.to_csv("path/to/GSE83654_correlation_matrix.csv")  # Save correlation matrix
    else:
        print("No LogRatio columns found to compute correlation.")


    ### Filtering for Significantly Regulated Genes ###
    threshold = 0.09 #Set a threshold
    num_datasets = len(logratio_cols)
    # Create boolean matrices for upregulation and downregulation
    up_matrix = pd.DataFrame({col: (merged_data[col] > threshold) for col in logratio_cols})
    down_matrix = pd.DataFrame({col: (merged_data[col] < -threshold) for col in logratio_cols})
    # Count how many datasets each gene is upregulated or downregulated
    up_mask = up_matrix.sum(axis=1) >= num_datasets
    down_mask = down_matrix.sum(axis=1) >= num_datasets
    combined_mask = up_mask | down_mask

    filtered_genes = merged_data[combined_mask]
    filtered_genes.to_csv(output_filtered, index=False)
    print("Filtered genes saved to GSE83654_filtered_genes.csv")

 ### Filtering Genes consistently upregulated or downregulated over time in all 3 cell lines ###
    print("\nFinding genes consistently upregulated or downregulated over time (PC3, DU145, LNCAP)")
    # Define the time-points columns for each cell line
    pc3_cols = ['_PC3_8h', '_PC3_24h', '_PC3_48h', '_PC3_72h']
    du145_cols = ['_DU145_8h', '_DU145_24h', '_DU145_48h', '_DU145_72h']
    lncap_cols = ['_LNCAP_8h', '_LNCAP_24h', '_LNCAP_48h', '_LNCAP_72h']

    # Filter for genes consistently upregulated over time in all three cell lines
    upregulated_trend = merged_data[
        (merged_data[pc3_cols[0]] < merged_data[pc3_cols[1]]) &
        (merged_data[pc3_cols[1]] < merged_data[pc3_cols[2]]) &
        (merged_data[pc3_cols[2]] < merged_data[pc3_cols[3]]) &

        (merged_data[du145_cols[0]] < merged_data[du145_cols[1]]) &
        (merged_data[du145_cols[1]] < merged_data[du145_cols[2]]) &
        (merged_data[du145_cols[2]] < merged_data[du145_cols[3]]) &

        (merged_data[lncap_cols[0]] < merged_data[lncap_cols[1]]) &
        (merged_data[lncap_cols[1]] < merged_data[lncap_cols[2]]) &
        (merged_data[lncap_cols[2]] < merged_data[lncap_cols[3]])
        ]

    # For downregulated over time
    downregulated_trend = merged_data[
        (merged_data[pc3_cols[0]] > merged_data[pc3_cols[1]]) &
        (merged_data[pc3_cols[1]] > merged_data[pc3_cols[2]]) &
        (merged_data[pc3_cols[2]] > merged_data[pc3_cols[3]]) &

        (merged_data[du145_cols[0]] > merged_data[du145_cols[1]]) &
        (merged_data[du145_cols[1]] > merged_data[du145_cols[2]]) &
        (merged_data[du145_cols[2]] > merged_data[du145_cols[3]]) &

        (merged_data[lncap_cols[0]] > merged_data[lncap_cols[1]]) &
        (merged_data[lncap_cols[1]] > merged_data[lncap_cols[2]]) &
        (merged_data[lncap_cols[2]] > merged_data[lncap_cols[3]])
        ]

    upregulated_trend.to_csv(output_up, index=False)
    downregulated_trend.to_csv(output_down, index=False)

    print(f"\nThe genes upregulated over time are saved to: {output_up} ({len(upregulated_trend)} genes)")
    print(f"The genes downregulated over time are saved to: {output_down} ({len(downregulated_trend)} genes)")
