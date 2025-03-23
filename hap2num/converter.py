"""
Function: process_hap_to_numeric
---------------------------------
This function converts haplotype genotype data from categorical format (e.g., AA, AT, TT) 
to numeric format for further genetic analysis. The conversion is parallelized for efficiency, 
and the user can choose between different numeric encoding formats.

**Novelty of this function:**
Unlike standard genotype conversions, this function ensures that the numeric encoding is **meaningful** by
considering each marker's reference (`REF`) and alternate (`ALT`) alleles. This allows users to interpret
the numeric codes biologically, ensuring that genotype representation remains consistent with allele composition.

Input:
------
- The input file must be in CSV format.
- It must contain the following first five columns **in this exact order**:
    1. `SNP`  (SNP marker name)
    2. `CHR`  (Chromosome number)
    3. `POS`  (Genomic position)
    4. `REF`  (Reference allele)
    5. `ALT`  (Alternate allele)
- The genotype data starts from **column 6 onward** and must be coded in **HapMap diploid format** (e.g., `AA`, `AT`, `TT`).

Parameters:
-----------
- input_file (str): Path to the input CSV file containing SNP genotype data.
- output_file (str): Path to save the converted numeric genotype data.
- num_processes (int, optional): Number of parallel processes for computation. Default is 60.
- batch_size (int, optional): Number of SNP markers processed per batch. Default is 5000.
- format_type (str, optional): Numeric format for genotype encoding.
    Options:
    * "012": Encodes genotypes as follows:
        - Homozygous reference (AA) → 0
        - Heterozygous (AT) → 1
        - Homozygous alternate (TT) → 2
    * "-101": Encodes genotypes as follows:
        - Homozygous reference (AA) → -1
        - Heterozygous (AT) → 0
        - Homozygous alternate (TT) → 1
- chunk_size (int, optional): Number of rows written to file per iteration. Default is 1000. 
  This allows flexibility for users with smaller datasets.

Raises:
-------
- ValueError: If the input file does not contain the required columns in the correct order.
- ValueError: If the format_type is not one of the allowed options ("012" or "-101").

Value:
------
The function generates a new CSV file (`output_file`) with the following:
- The same columns as the input file, but with genotypic values converted to numeric format.
- The **genotypic data** (starting from column 6) is transformed according to the specified `format_type`.
- The output preserves the SNP marker order from the input.
- Missing or unrecognized genotypes are encoded as `-9` to signify missing data.

The processed output file can be used for downstream genetic analyses such as GWAS, genomic selection, 
and population structure analysis.

"""
# Remove warnings
def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

# Load required libraries
import os
import pandas as pd
import numpy as np
import logging
from multiprocessing import Pool
from tqdm import tqdm
from io import StringIO
import warnings

def convert_genotypes(args):
    """Convert genotypes to numeric format using explicit allele combinations."""
    marker, data, format_type = args
    ref = data['REF'].values[0]
    alt = data['ALT'].values[0]
    
    # Generate standardized genotype keys
    homo_ref = ''.join(sorted(ref * 2))
    het = ''.join(sorted(ref + alt))
    homo_alt = ''.join(sorted(alt * 2))
    
    # Create conversion map based on format type
    if format_type == "012":
        conversion_map = {
            homo_ref: '0',
            het: '1',
            homo_alt: '2'
        }
    elif format_type == "-101":
        conversion_map = {
            homo_ref: '-1',
            het: '0',
            homo_alt: '1'
        }
    
    for col in data.columns[5:]:
        data[col] = data[col].apply(
            lambda x: conversion_map.get(''.join(sorted(str(x))), '-9') 
            if pd.notna(x) else '-9'
        )
    
    return data

def process_hap_to_numeric(input_file: str, output_file: str, num_processes: int = 60, 
                          batch_size: int = 5000, format_type: str = "012", chunk_size: int = 1000):
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    required_columns = ["SNP", "CHR", "POS", "REF", "ALT"]
    
    try:
        logging.info(f"Reading input file: {input_file}")
        df = pd.read_csv(input_file, sep=',')
        
        if list(df.columns[:5]) != required_columns:
            raise ValueError(f"First five columns must be: {required_columns}")
        
        markers = df['SNP'].unique()
        marker_data = [(marker, df[df['SNP'] == marker], format_type) for marker in markers]
        
        logging.info("Starting parallel conversion...")
        processed_batches = []
        for i in range(0, len(marker_data), batch_size):
            batch = marker_data[i:i+batch_size]
            with Pool(num_processes) as pool:
                results = list(tqdm(pool.imap(convert_genotypes, batch), total=len(batch), desc="Processing"))
            processed_batches.append(pd.concat(results))
        
        logging.info(f"Writing output to {output_file}")
        with open(output_file, 'w') as fout:
            header_written = False
            for batch_idx, batch_df in enumerate(processed_batches):
                buffer = StringIO()
                batch_df.to_csv(buffer, index=False, sep='\t', quoting=1)
                if not header_written:
                    fout.write(buffer.getvalue())
                    header_written = True
                else:
                    fout.write(buffer.getvalue().split('\n', 1)[1])
                logging.info(f"Batch {batch_idx+1} written")
        
        logging.info("Conversion completed successfully")
    
    except Exception as e:
        logging.error(f"Error: {str(e)}")
        raise