# Import
import pandas as pd
import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import PandasTools,Draw
import math
import openbabel
import darkchem
import sys
sys.path.append("../..")
import functions

# Define testing function
def test_struc2mol():
    # Generate some data to use
    sms = {'smiles': ['CC(=O)C', 'O=C1CCCCC1', ' ', 'O=C1CCCC1C', 'CCCC(=O)CCC',
                          'CCCCC(=O)CCCC', 'CC(=O)CC(=O)C', ' ', 'COCC(=O)C',
                          'COCC(=O)C'] }
    sms = pd.DataFrame(sms)

    # Operate the function wihch is being tested
    struc2mol = functions.struc2mol(sms)

    # Assert several arguments
    assert len(struc2mol.loc[0]) == 3
    assert type(struc2mol.loc[2][2]) == type(struc2mol.loc[7][2])

    return 0
