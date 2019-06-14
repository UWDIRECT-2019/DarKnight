import pandas as pd
import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import PandasTools,Draw
import math
import openbabel
import darkchem
from .. import fcts

def test_struc2mol():
    # Generate some data to use
    sms = {'smiles': ['CC(=O)C', 'O=C1CCCCC1', ' ', 'O=C1CCCC1C', 'CCCC(=O)CCC',
                          'CCCCC(=O)CCCC', 'CC(=O)CC(=O)C', ' ', 'COCC(=O)C',
                          'COCC(=O)C'] }
    sms = pd.DataFrame(sms)

    struc2mol = fcts.struc2mol(sms)

    assert len(struc2mol.loc[0]) == 3
    assert type(struc2mol.loc[2][2]) == type(struc2mol.loc[7][2])

    return 0
