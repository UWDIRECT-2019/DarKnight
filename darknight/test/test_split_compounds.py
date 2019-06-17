import pandas as pd
import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import PandasTools,Draw
import math
import openbabel
import darkchem
from IPython.display import display
from .. import fcts

def test_split_compounds():
    # Import some data to use
    data = pd.read_excel('test_split_compounds.xlsx')

    split_compounds = fcts.split_compounds(data)

    assert len(data) == len(split_compounds)
    assert type(data) == type(split_compounds)
    assert split_compounds.shape == (5,2)

    return 0
