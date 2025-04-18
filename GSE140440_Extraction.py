# ===================================================================================== #
#  WRITTEN BY : Asmaa Yassin                                                            #
#  DATE       : April 17, 2025                                                          #
# ------------------------------------------------------------------------------------- #
#  DESCRIPTION: Extracts all .gz files in a directory and saves the uncompressed files. #
# ===================================================================================== #


import gzip
import shutil
import os

directory = "path/to/GSE######_RAW"  # Set the correct directory where .gz count files are located

for filename in os.listdir(directory): # Loop through all .gz files and extract them
    if filename.endswith(".gz"):
        filepath = os.path.join(directory, filename)
        output_filepath = filepath[:-3]  # Removes ".gz" extension

        with gzip.open(filepath, 'rb') as f_in: #decompress the files and save the uncompressed to a new file
            with open(output_filepath, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        print(f"Extracted: {filename} -> {output_filepath}") #print a message showing the input and output file names
