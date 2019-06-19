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
def test_remove_other_elements():
    # Import some data to use
    data = pd.read_excel('test_remove_other_elements.xlsx')

    # Operate the function wihch is being tested
    remove_other_elements = process.remove_other_elements(data)

    # Assert several arguments
    assert len(remove_other_elements.loc[remove_other_elements['Reactants'] == 'F']) == len(remove_other_elements.loc[remove_other_elements['Reactants'] == 'l'])
    assert len(remove_other_elements.loc[remove_other_elements['Reactants'] == 'B']) == len(remove_other_elements.loc[remove_other_elements['Reactants'] == 'r'])
    assert len(remove_other_elements.loc[remove_other_elements['Reactants'] == 'I']) == len(remove_other_elements.loc[remove_other_elements['Reactants'] == 'i'])
    assert len(remove_other_elements.loc[remove_other_elements['Reactants'] == 'M']) == len(remove_other_elements.loc[remove_other_elements['Reactants'] == 'g'])
    assert len(remove_other_elements.loc[remove_other_elements['Reactants'] == 'L']) == len(remove_other_elements.loc[remove_other_elements['Reactants'] == 'b'])
    assert len(remove_other_elements.loc[remove_other_elements['Reactants'] == 'a']) == len(remove_other_elements.loc[remove_other_elements['Reactants'] == 'e'])
    assert len(remove_other_elements.loc[remove_other_elements['Reactants'] == 'K']) == len(remove_other_elements.loc[remove_other_elements['Reactants'] == 'V'])
    assert len(remove_other_elements.loc[remove_other_elements['Reactants'] == 'd']) == len(remove_other_elements.loc[remove_other_elements['Reactants'] == 'R'])
    assert len(remove_other_elements.loc[remove_other_elements['Reactants'] == 'Z']) == len(remove_other_elements.loc[remove_other_elements['Reactants'] == 'G'])
    assert len(remove_other_elements.loc[remove_other_elements['Reactants'] == 'A']) == len(remove_other_elements.loc[remove_other_elements['Reactants'] == 'Y'])
    assert len(remove_other_elements.loc[remove_other_elements['Reactants'] == 'u']) == 0

    return 0
