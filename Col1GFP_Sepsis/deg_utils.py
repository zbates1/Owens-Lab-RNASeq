import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np


def read_differential_expression_files(deg_filename):
    """
    Reads an Excel file and returns a Pandas DataFrame
    
    Parameters
    ----------
    filename : str
        Path to Excel file
    
    Returns
    -------
    df : pd.DataFrame
        DataFrame containing data from Excel file
    """
    # assert file is of excel type
    assert deg_filename.endswith('.xlsx'), f'File {deg_filename} is not an Excel file'
    # assert 
    df = pd.read_excel(deg_filename)
    print(df.head())

    return df

def read_mrna_gene_expression_files(mrna_filename):
    """
    Reads an Excel file containing mRNA gene expression data (weighted average based on transcript length already computed)and returns a Pandas DataFrame
    
    Parameters
    ----------
    filename : str
        Path to Excel file
    
    Returns
    -------
    df : pd.DataFrame
        DataFrame containing data from Excel file
    """

    if mrna_filename.endswith('.xls') or mrna_filename.endswith('.xlsx'):
        df = pd.read_excel(mrna_filename)  # Assuming it's a real Excel file
    elif mrna_filename.endswith('.csv'):
        df = pd.read_csv(mrna_filename)  # For CSV files
    elif mrna_filename.endswith('.xls') and not mrna_filename.endswith('.xlsx'):
        # Attempt to read as a TSV if .xls files are actually mislabeled
        try:
            df = pd.read_csv(mrna_filename, sep='\t')  # Try reading as a tab-separated file
        except Exception as e:
            raise ValueError(f"Failed to read the file with error: {e}")
    else:
        raise ValueError("Unsupported file extension. Please provide a .csv, .xls, or .xlsx file.")
    print(df.head())

    return df

def compute_transcript_avg_counts_conditions(filename):
    """
    Computes transcript average counts for each condition:
    - D5
    - Control
    - D14
    

    The Mice are grouped in the following conditions:
    Sample ID,Condition,Group
    237N,Control,14d
    281N,Control,14d
    267N,Sepsis,5d
    284L,Sepsis,5d
    284N,Sepsis,5d
    280N,Sepsis,5d
    237R,Sepsis,14d
    238R,Sepsis,14d
    266R,Sepsis,14d


    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing data from Excel file
    
    Returns
    -------
    df : pd.DataFrame
        DataFrame containing data from Excel file
    """

    # Gather the df from the xls file
    df = read_mrna_gene_expression_files(filename)


    # Tag: ASK OWEN -> one fo the D5 samples is a control in the mRNA_gene.FPKM.csv file??
    # d14_cohort = ['237R', '238R', '266R']
    # d5_cohort = ['267N', '284L', '284N', '280N']
    # cntl_cohort = ['237N', '281N']

    # New sample names based on the csv file
    d14_cohort = ['D14_1_237R', 'D14_2_238R', 'D14_3_266R']
    d5_cohort = ['D5_1_267N', 'D5_2_284L', 'D5_3_284N']
    cntl_cohort = ['C_1_237N', 'C_2_281N', 'C_3_280N']

    # Initialize dictionaries to store DataFrame slices
    d14_mean_df = {}
    d5_mean_df = {}
    cntl_mean_df = {}

    # Extract DataFrames for each cohort
    for sample in d14_cohort:
        d14_mean_df[sample] = df[sample]
    for sample in d5_cohort:
        d5_mean_df[sample] = df[sample]
    for sample in cntl_cohort:
        cntl_mean_df[sample] = df[sample]

    # Convert dictionaries to DataFrames and calculate the mean across the columns
    d14_df = pd.DataFrame(d14_mean_df).mean(axis=1)
    d5_df = pd.DataFrame(d5_mean_df).mean(axis=1)
    cntl_df = pd.DataFrame(cntl_mean_df).mean(axis=1)

    # Merge the mean DataFrames into a final DataFrame
    final_df = pd.DataFrame({
        'D14': d14_df,
        'D5': d5_df,
        'Control': cntl_df
    })

    # Add the Gene Name column if it exists in the original df
    if 'gene_name' in df.columns:
        final_df['Gene Name'] = df['gene_name']

    print(f'\nFinal df: {final_df.head()}\n')
    print(f'Final df shape: {final_df.shape}')

    return final_df
