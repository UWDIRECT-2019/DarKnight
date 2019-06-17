import pandas as pd
import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import PandasTools,Draw
import math
import openbabel
import darkchem
from IPython.display import display
from .. import fcts

def test_data_clean():
    # Import some data to use
    data = pd.read_excel('test_data_clean.xlsx')

    data_clean = fcts.data_clean(data)

    assert len(data_clean.loc[data_clean['Reactants'] == 'F']) == len(data_clean.loc[data_clean['Reactants'] == 'l'])
    assert len(data_clean.loc[data_clean['Reactants'] == 'B']) == len(data_clean.loc[data_clean['Reactants'] == 'r'])
    assert len(data_clean.loc[data_clean['Reactants'] == 'I']) == len(data_clean.loc[data_clean['Reactants'] == 'i'])
    assert len(data_clean.loc[data_clean['Reactants'] == 'M']) == len(data_clean.loc[data_clean['Reactants'] == 'g'])
    assert len(data_clean.loc[data_clean['Reactants'] == 'L']) == len(data_clean.loc[data_clean['Reactants'] == 'b'])
    assert len(data_clean.loc[data_clean['Reactants'] == 'a']) == len(data_clean.loc[data_clean['Reactants'] == 'e'])
    assert len(data_clean.loc[data_clean['Reactants'] == 'K']) == len(data_clean.loc[data_clean['Reactants'] == 'V'])
    assert len(data_clean.loc[data_clean['Reactants'] == 'd']) == len(data_clean.loc[data_clean['Reactants'] == 'R'])
    assert len(data_clean.loc[data_clean['Reactants'] == 'Z']) == len(data_clean.loc[data_clean['Reactants'] == 'G'])
    assert len(data_clean.loc[data_clean['Reactants'] == 'A']) == len(data_clean.loc[data_clean['Reactants'] == 'Y'])
    assert len(data_clean.loc[data_clean['Reactants'] == 'u']) == 0

    return 0
