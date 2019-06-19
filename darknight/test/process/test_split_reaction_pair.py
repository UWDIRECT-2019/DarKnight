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
import process

# Define testing function
def test_split_reaction_pair():
    # Import some data to use
    data = pd.read_excel('test_split_reaction_pair.xlsx')

    # Operate the function wihch is being tested
    split_reaction_pair = process.split_reaction_pair(data)

    # Assert several arguments
    assert len(data) == len(split_reaction_pair)
    assert type(data) == type(split_reaction_pair)
    assert split_reaction_pair.shape == (5,2)

    return 0
