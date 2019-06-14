import pandas as pd
import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import PandasTools,Draw
import math
import openbabel
import darkchem
from .. import fcts

def test_vector_magnitude():
    # Load model
    model = darkchem.utils.load_model('/Family/David/UW_ChemE/2019_Spring/ChemE_547/DarKnight_File/Final Trained DarkChem Network Weights/N7b_[M+H]/')
    # Generate some data to use
    smilescc = {'reactants': ['CC(=O)C', 'O=C1CCCCC1', 'O=C1CCCC1C', 'O=C1CCCC1C', 'CCCC(=O)CCC',
                          'CCCCC(=O)CCCC', 'CC(=O)CC(=O)C', 'CC(=O)CC(=O)C', 'COCC(=O)C',
                          'COCC(=O)C',
                          'CC(=O)c1ccccc1']}
    # Conversion to DarkChem latent space vectors
    smilescc['rvec'] = [darkchem.utils.struct2vec(reactant).astype(np.int16) for reactant in smilescc['reactants']]
    smilescc['rlat'] = model.encoder.predict(np.array(smilescc['rvec']))
    # Coversion to DataFrame
    smilescc['rlat'] = pd.DataFrame(smilescc['rlat'])

    vector_magnitude = fcts.vector_magnitude(smilescc['rlat'])

    NoneType = type(None)
    x = None

    assert isinstance(vector_magnitude, NoneType)

    return 0
