import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np


def read_excel(filename):
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
    assert filename.endswith('.xlsx'), f'File {filename} is not an Excel file'
    # assert 
    df = pd.read_excel(filename)
    print(df.head())