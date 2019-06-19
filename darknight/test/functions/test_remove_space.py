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
def test_remove_space():
    # Generate some data to use
    data = {'Reactants': ['C /C=C/C(=O)C#N', 'C CNCC(=O)C'],
             'Products': ['C /C=C/C(C#N)O', 'C CNCC(O)C']}
    data1 = {'Reactants': ['C/C=C/C(=O)C#N', 'CCNCC(=O)C'],
             'Products': ['C/C=C/C(C#N)O', 'CCNCC(O)C']}
    data = pd.DataFrame(data)
    data1 = pd.DataFrame(data1)

    # Operate the function wihch is being tested
    data = functions.remove_space(data)

    # Assert several arguments
    assert len(data['Reactants'][0]) == len(data1['Reactants'][0])
    assert len(data['Reactants'][1]) == len(data1['Reactants'][1])
    assert len(data['Products'][0]) == len(data1['Products'][0])
    assert len(data['Products'][1]) == len(data1['Products'][1])

    return 0
