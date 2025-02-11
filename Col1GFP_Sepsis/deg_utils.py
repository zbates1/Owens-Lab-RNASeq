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
    Reads an Excel file containing mRNA gene expression data and returns a Pandas DataFrame.
    """
    if mrna_filename.endswith(('.xls', '.xlsx')):
        df = pd.read_excel(mrna_filename)
    elif mrna_filename.endswith('.csv'):
        df = pd.read_csv(mrna_filename)
    else:
        raise ValueError("Unsupported file extension. Please provide a .csv, .xls, or .xlsx file.")

    # Ensure consistent column naming
    if 'gene_name' in df.columns:
        df.rename(columns={'gene_name': 'Gene Name'}, inplace=True)

    print(f'Shape of df: {df.shape}')
    print(f'DF Columns after rename: {df.columns}')
    print(df.head())
    return df

def compute_transcript_avg_counts_conditions(filename):
    """
    Computes transcript average counts for each condition: D5, Control, D14.
    """
    df = read_mrna_gene_expression_files(filename)

    # Define cohorts based on sample naming in your data
    d14_cohort = ['D14_1_237R', 'D14_2_238R', 'D14_3_266R']
    d5_cohort = ['D5_1_267N', 'D5_2_284L', 'D5_3_284N']
    cntl_cohort = ['C_1_237N', 'C_2_281N', 'C_3_280N']

    # Use pandas to directly calculate the mean of the cohorts
    d14_mean = df[d14_cohort].mean(axis=1)
    d5_mean = df[d5_cohort].mean(axis=1)
    cntl_mean = df[cntl_cohort].mean(axis=1)

    # Creating a new DataFrame from the means
    final_df = pd.DataFrame({
        'D14': d14_mean,
        'D5': d5_mean,
        'Control': cntl_mean
    })

    # Assuming 'Gene Name' or some identifier is in the original df to use as an index
    final_df['Gene Name'] = df['Gene Name']

    # Set 'Gene Name' as index if it's intended to be the identifier

    print(f'\nFinal df: {final_df.head()}\n')
    print(f'Final df shape: {final_df.shape}')

    return final_df
