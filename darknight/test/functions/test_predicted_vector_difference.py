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

def test_predicted_vector_difference():
    # Load model
    model = darkchem.utils.load_model('/Family/David/UW_ChemE/2019_Spring/ChemE_547/DarKnight_File/Final Trained DarkChem Network Weights/N7b_[M+H]/')
    # Generate some data to use
    smilescc = {'Reactants': ['CC(=O)C', 'O=C1CCCCC1', 'O=C1CCCC1C', 'O=C1CCCC1C', 'CCCC(=O)CCC',
                          'CCCCC(=O)CCCC', 'CC(=O)CC(=O)C', 'CC(=O)CC(=O)C', 'COCC(=O)C',
                          'COCC(=O)C',
                          'CC(=O)c1ccccc1'],
             'Products': ['CCC(=O)C', 'CC1CCCCC1=O', 'CC1CCC(C1=O)C', 'O=C1CCCC1(C)C', 'CCCC(=O)C(CC)C',
                          'CCCCC(=O)C(CCC)C', 'CC(C(=O)C)C(=O)C', 'CC(=O)C(C(=O)C)(C)C', 'COCC(=O)CC',
                          'COC(C(=O)CC)C',
                          'CCC(=O)c1ccccc1']}
    # Conversion to DarkChem latent space vectors
    smilescc['rvec'] = [darkchem.utils.struct2vec(reactant).astype(np.int16) for reactant in smilescc['Reactants']]
    smilescc['pvec'] = [darkchem.utils.struct2vec(product).astype(np.int16) for product in smilescc['Products']]
    smilescc['rlat'] = model.encoder.predict(np.array(smilescc['rvec']))
    smilescc['plat'] = model.encoder.predict(np.array(smilescc['pvec']))
    # Coversion to DataFrame
    smilescc['rlat'] = pd.DataFrame(smilescc['rlat'])
    smilescc['plat'] = pd.DataFrame(smilescc['plat'])

    # Operate the function wihch is being tested
    difference = functions.predicted_vector_difference(smilescc['rlat'],smilescc['plat'])

    # Assert several arguments
    assert isinstance(difference, list)
    assert len(difference) == 11

    return 0
