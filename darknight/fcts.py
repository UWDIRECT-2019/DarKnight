import pandas as pd
import numpy as np
from rdkit.Chem import AllChem as Chem


def array_in_nd_array(test, array):
    """
    Checks whether or not a test 1D array is contained within a full ND array.
    Returns True if the test array is equal to any of the dimensions of the ND array.
    Returns False if the test array does not match any dimension of the ND array.
    """
    return any(np.array_equal(x, test) for x in array)

# A function to remove the space in the string
def remove_space(data):
    """
    Remove the intermediate redundant space in the smiles strings,
    The name of the column must be 'Reactants' and 'Products'
    
    """
    for i in range(data.shape[0]):
        data['Reactants'][i] = data['Reactants'][i].replace(' ','')
        data['Products'][i] = data['Products'][i].replace(' ','')
    return data

# the function to count Pearson correlation
def r2pcorr(data1,data2):
    """
    A function to calculate the Pearson Correlation Coefficient
    between the latent space vectors of reactants and products.
    """
    metric = pd.DataFrame(columns = ['Correlation'])
    for i in range(data1.shape[0]):
        metric.loc[i,'Correlation'] = data1.iloc[i].corr(data2.iloc[i])
    return metric

def struc2mol(sms):
    """
    A function to transform smiles strings to molecules with the module 
    rdkit.Chem.MolFromSmiles, and return a DataFrame
    """
    save = pd.DataFrame(columns = ['smiles','mol'])
    save['smiles'] = sms['smiles']
    for i in range(sms.shape[0]):
        save['mol'][i] = Chem.MolFromSmiles(sms['smiles'][i])
    return save

def data_clean(data):
    """
    screen chemicla reactions only possess C/H/O/N/P/S elements
    """
    charset = ['F','l','B','r','I','i','M','g','L','b','a','e','K','V','d','R','Z','G','A','Y','u','H']
    x = []
    for i in range(data.shape[0]):
        for j in range(len(data.iloc[i,1])):
            if data.iloc[i,1][j] in charset:
                x.append(i)
                break
    df = data[(True^data['Index'].isin(x))]
    df.reset_index(drop=True, inplace=True)
    return df

def split_compounds(data):
    """
    split the reactants and prodcuts within one complete chemical reaction
    """
    df = pd.DataFrame(columns = ['Reactants','Products'])
    for i in range(data.shape[0]):
        a = data.iloc[i][1].split('	')
        df.loc[i,'Reactants'] = a[0]
        df.loc[i,'Products'] = a[1]
    return df