"""
**hap2num**

`hap2num` is a Python package for converting haplotype genotype data into numeric format. It ensures meaningful numerical representation based on allele composition.

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

## Installation

To install from GitHub:
```bash
pip install git+https://github.com/FAkohoue/hap2num.git
