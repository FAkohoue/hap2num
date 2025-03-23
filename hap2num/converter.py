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
- format_type (str, optional): Numeric format for genotype encoding.
    Options:
    * "012": Encodes genotypes as follows:
        - Homozygous reference (AA) → 0
        - Heterozygous (AT) → 1
        - Homozygous alternate (TT) → 2
    * "-101": Encodes genotypes as follows:
        - Homozygous reference (AA) → 1
        - Heterozygous (AT) → 0
        - Homozygous alternate (TT) → -1
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
# Load required libraries
import os
import pandas as pd
import numpy as np
import logging
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
from functools import partial
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings("ignore")

# Create conversion map based on format type
def get_binary_genotype_map(ref_allele, alt_allele, format_type: str = "012"):
    if format_type not in {"012", "-101"}:
        raise ValueError("Invalid format_type. Choose '012' or '-101'.")

    return {
        '012': {
            f'{ref_allele}{ref_allele}': '0',
            f'{ref_allele}{alt_allele}': '1',
            f'{alt_allele}{ref_allele}': '1',
            f'{alt_allele}{alt_allele}': '2'
        },
        '-101': {
            f'{ref_allele}{ref_allele}': '1',
            f'{ref_allele}{alt_allele}': '0',
            f'{alt_allele}{ref_allele}': '0',
            f'{alt_allele}{alt_allele}': '-1'
        }
    }[format_type]

def process_row(row, format_type):
    """Process a single row with vectorized operations"""
    ref = row['REF']
    alt = row['ALT']
    genotype_map = get_binary_genotype_map(ref, alt, format_type)
    
    # Vectorized replacement for all sample columns
    sample_data = row[5:].replace(genotype_map)
    # Handle missing/unmapped values
    sample_data = sample_data.where(sample_data.isin(genotype_map.values()), '-9')
    
    return pd.concat([row[:5], sample_data])

def process_chunk(chunk, format_type):
    """Process a chunk of rows with vectorized row operations"""
    return chunk.apply(process_row, axis=1, format_type=format_type)

def process_hap_to_numeric(input_file: str, output_file: str, 
                          num_processes: int = None, 
                          chunk_size: int = 1000,
                          format_type: str = "012"):
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    # Validate format type before proceeding
    if format_type not in {"012", "-101"}:
        raise ValueError("Invalid format_type. Choose '012' or '-101'.")

    required_columns = ["SNP", "CHR", "POS", "REF", "ALT"]
    
    try:
        # Auto-detect CPU cores if not specified
        num_processes = num_processes or cpu_count()
        
        logging.info(f"Reading input file: {input_file}")
        df = pd.read_csv(input_file, sep=',', dtype={'CHR': str})
        
        if list(df.columns[:5]) != required_columns:
            raise ValueError(f"First five columns must be: {required_columns}")
        
        # Pre-convert sample columns to categorical type
        sample_cols = df.columns[5:]
        df[sample_cols] = df[sample_cols].astype('category')
        
        # Split data into row chunks for parallel processing
        total_rows = len(df)
        chunks = [df.iloc[i:i+chunk_size] for i in range(0, total_rows, chunk_size)]
        
        logging.info(f"Processing {total_rows} SNPs in {len(chunks)} chunks using {num_processes} cores")
        
        # Create partial function for static parameters
        process_chunk_partial = partial(process_chunk, format_type=format_type)
        
        with Pool(num_processes) as pool:
            results = []
            with tqdm(total=len(chunks), desc="Processing chunks") as pbar:
                for chunk_result in pool.imap(process_chunk_partial, chunks):
                    results.append(chunk_result)
                    pbar.update()
        
        logging.info("Combining results and optimizing memory")
        final_df = pd.concat(results, ignore_index=True)
        
        # Convert categorical back to string for output
        final_df[sample_cols] = final_df[sample_cols].astype(str)
        
        logging.info(f"Writing output to {output_file}")
        final_df.to_csv(output_file, index=False, encoding='utf-8')
        
        logging.info("Conversion completed successfully")
        return True
    
    except Exception as e:
        logging.error(f"Critical error: {str(e)}")
        raise RuntimeError(f"Processing failed: {str(e)}") from e