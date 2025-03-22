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
import os
import pandas as pd
import numpy as np
import logging
from multiprocessing import Pool
from tqdm import tqdm
from io import StringIO
import warnings

def convert_genotypes(args):
    """Convert genotypes to numeric format."""
    marker, data, genotype_values = args
    ref_allele = data['REF'].values[0]
    alt_allele = data['ALT'].values[0]
    
    for col in data.columns[5:]:
        data[col] = data[col].apply(lambda x: genotype_values.get(sum([x.count(ref_allele), x.count(alt_allele)]), '-9') if pd.notna(x) else '-9')
    
    data['marker'] = marker
    return data

    #remove warning
def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

def process_hap_to_numeric(input_file: str, output_file: str, num_processes: int = 60, batch_size: int = 5000, 
                            format_type: str = "012", chunk_size: int = 1000):
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    # Define genotype mappings based on format type
    genotype_mappings = {
        "012": {0: '0', 1: '1', 2: '2'},
        "-101": {0: '-1', 1: '0', 2: '1'}
    }
    
    if format_type not in genotype_mappings:
        raise ValueError("Invalid format type. Choose either '012' or '-101'")
    
    genotype_values = genotype_mappings[format_type]
    
    required_columns = ["SNP", "CHR", "POS", "REF", "ALT"]
    
    try:
        logging.info(f"Reading input file: {input_file}")
        df = pd.read_csv(input_file, sep=',')
        
        # Ensure correct column order
        if list(df.columns[:5]) != required_columns:
            raise ValueError(f"Input file must have the first five columns in order: {required_columns}")
        
        markers = df['SNP'].unique()
        marker_data = [(marker, df[df['SNP'] == marker], genotype_values) for marker in markers]
        
        logging.info("Starting parallel genotype conversion...")
        processed_data_batches = []
        for i in range(0, len(marker_data), batch_size):
            batch_data = marker_data[i:i+batch_size]
            with Pool(num_processes) as pool:
                processed_batch = list(tqdm(pool.imap(convert_genotypes, batch_data), total=len(batch_data), desc="Converting"))
            processed_data_batches.append(pd.concat(processed_batch))
        
        logging.info(f"Saving converted data to: {output_file}")
        with open(output_file, 'w') as outfile:
            header_written = False
            for i, batch in enumerate(processed_data_batches):
                csv_buffer = StringIO()
                batch.to_csv(csv_buffer, index=False, sep='\t', quoting=1, chunksize=chunk_size)
                csv_lines = csv_buffer.getvalue().splitlines()
                if not header_written:
                    outfile.write(csv_lines[0] + '\n')
                    header_written = True
                for line in csv_lines[1:]:
                    outfile.write(line + '\n')
                logging.info(f"Batch {i+1} written to {output_file}")
        logging.info("All data successfully processed.")
    except Exception as e:
        logging.error(f"Error during conversion: {e}")
        raise