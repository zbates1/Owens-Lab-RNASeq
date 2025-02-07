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
    # assert file is of excel type
    assert mrna_filename.endswith('.xls'), f'File {mrna_filename} is not an xls file'
    # assert 
    df = pd.read_excel(mrna_filename)
    print(df.head())

    return df

def compute_transcript_avg_counts_conditions(df, filename):
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

    d14_cohort = ['237R', '238R', '266R']
    d5_cohort = ['267N', '284L', '284N', '280N']
    cntl_cohort = ['237N', '281N']
    