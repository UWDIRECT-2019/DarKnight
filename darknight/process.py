"""
darknight.process
~~~~~~~~~~~~~~~~~
Contains most of the data cleaning functions built for Darknight.
"""

# Imports

import pandas as pd
import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import PandasTools,Draw
import math
import openbabel
import darkchem


# Functions

def remove_other_elements(data):
    """Screens for chemical reactions that only contain C/H/O/N/P/S elements
    """
    charset = ['F','l','B','r','I','i','M','g','L','b','a','e','K','V','d','R','Z','G','A','Y','u']
    x = []
    for i in range(data.shape[0]):
        for j in range(len(data.iloc[i,1])):
            if data.iloc[i,1][j] in charset:
                x.append(i)
                break
    df = data[(True^data['Index'].isin(x))]
    df.reset_index(drop=True, inplace=True)
    return df

def split_reaction_pair(data):
    """Takes a dataframe with Reactants and Products in separate columns.
    The Reactants and Products represent a chemical reaction pair.
    Splits the reactants and products within one complete chemical reaction.
    """
    df = pd.DataFrame(columns = ['Reactants','Products'])
    for i in range(data.shape[0]):
        a = data.iloc[i][1].split('→')
        df.loc[i,'Reactants'] = a[0]
        df.loc[i,'Products'] = a[1]
    return df
