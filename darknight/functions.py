"""
DarKnight.functions
~~~~~~~~~~~~~~~~~~~
General utility functions for DarKnight.
"""

# Imports

import pandas as pd
import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import PandasTools, Draw
import math
import openbabel
import darkchem
from IPython.display import display
import tensorflow as tf
tf.logging.set_verbosity(tf.logging.ERROR)


# Functions

def array_in_nd_array(test, array):
    """
    Checks whether or not a test 1D array is contained within a full ND array.
    Returns True if the test array is equal to any of the dimensions of the ND array.
    Returns False if the test array does not match any dimension of the ND array.
    """
    return any(np.array_equal(x, test) for x in array)



def remove_space(data):
    """Removes the intermediate redundant space in SMILES strings.
    The input must contain two columns titled 'Reactants' and 'Products'

    """
    for i in range(data.shape[0]):
        data['Reactants'][i] = data['Reactants'][i].replace(' ', '')
        data['Products'][i] = data['Products'][i].replace(' ', '')
    return data


def r2pcorr(data1, data2):
    """A function to calculate the Pearson correlation coefficient
    between the latent space vectors of reactants and products.
    """
    metric = pd.DataFrame(columns = ['Correlation'])
    for i in range(data1.shape[0]):
        metric.loc[i,'Correlation'] = data1.iloc[i].corr(data2.iloc[i])
    return metric


def struc2mol(sms):
    """A function to transform SMILES strings to molecules with the module
    rdkit.Chem.MolFromSmiles, and returns DataFrame
    """
    save = pd.DataFrame(columns = ['raw_smiles','smiles','mol'])
    save['raw_smiles'] = sms['smiles']
    for i in range(sms.shape[0]):
        save['mol'][i] = Chem.MolFromSmiles(sms['smiles'][i])
        if save['mol'][i] is None:
            save['smiles'][i] = 'Invalid smi str'
        else:
            save['smiles'][i] = sms['smiles'][i]
    return save


def predicted_vector_difference(actualvec, predvec):
    """Calculates the difference between actual and predicted
    product vectors in latent space.
    """
    d = []
    for i in range(len(actualvec)):
        s = 0
        for j in range(actualvec.shape[1]):
            s += (actualvec.iloc[i][j] - predvec.iloc[i][j])**2
        s = np.sqrt(s)
        d.append(s)
    return d


def vector_magnitude(data):
    """Computes the average and standard deviation of the reaction path vector
    magnitudes.
    """
    a = []
    for i in range(len(data)):
        s = 0
        for j in range(data.shape[1]):
            s += (data.iloc[i][j])**2
        s = np.sqrt(s)
        a.append(s)
    avg = np.average(a)
    std = np.std(a)
    print ('The average magnitude is:', avg)
    print ('The standard deviation is:', std)
    #return avg, std


def vector_angle(rct, prd):
    """Computes the average and standard deviation of the reaction path vector
    angles.
    """
    #u = []
    #d = []
    angle = []
    for i in range(len(rct)):
        up = 0
        rm = 0
        pm = 0
        for j in range(rct.shape[1]):
            up += rct.iloc[i][j] * prd.iloc[i][j]  #numerator
            rm += (rct.iloc[i][j])**2  # the magnitude of reactant vector
            pm += (prd.iloc[i][j])**2  # the magnitude of product vector
        #u.append(up)
        rm = np.sqrt(rm)
        pm = np.sqrt(pm)
        cos = up/(rm*pm)
        a = math.degrees(math.acos(cos))
        #d.append(rm*pm)
        angle.append(a)
    aveg = np.average(angle)
    std = np.std(angle)
    print('The average angle is:', aveg)
    print('The standard deviation is:', std)


def standardize_smi(smi):
    """Standardizes SMILES strings into SMILES strings through
    OpenBabel.
    (Important in optimizing prediction results.)
    """
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi", "smi")
    mol = openbabel.OBMol()
    obConversion.ReadString(mol, smi)
    outMDL = obConversion.WriteString(mol)[:-2]
    return outMDL

def standardize_can(smi):
    """Standardizes SMILES strings into canonical SMILES strings through
    OpenBabel.
    (Important in optimizing prediction results.)
    """
    obconversion = openbabel.OBConversion()
    obconversion.SetOutFormat("can")
    #obconversion.SetInAndOutFormats("smi", "can")
    obmol = openbabel.OBMol()
    obconversion.ReadString(obmol, smi)
    outMDL = obconversion.WriteString(obmol)[:-2]
    return outMDL


def path_vec(data, model):
    """Calculates the reaction path vector for each type of chemical reaction.

    Args:
        data: Dataframe with 'Reactants' and 'Products' as columns
        model: Trained model containing latent space

    Returns:
        A 128-dimensional numpy array representation of the mean path vector
            between reactants and products contained in the input dataframe
    """
    rvec = [darkchem.utils.struct2vec(reactant) for reactant in data['Reactants']]
    pvec = [darkchem.utils.struct2vec(product) for product in data['Products']]
    rvec = np.array(rvec).astype(int)
    pvec = np.array(pvec).astype(int)
    r_latent = model.encoder.predict(rvec)
    p_latent = model.encoder.predict(pvec)
    rvecdf = pd.DataFrame(r_latent)
    pvecdf = pd.DataFrame(p_latent)
    path = pvecdf - rvecdf
    path_vec = np.array(path.mean().values)
    return path_vec


def transform_r2p_str(smi, model, path_vec, k):
    """Transforms reactant SMILES string to product SMILES string
    """
    test = darkchem.utils.struct2vec(smi)
    test = np.array(test)
    test = test.reshape(-1,100)
    t_l = model.encoder.predict(test)
    t_pre = t_l + path_vec
    t_pred = model.decoder.predict(t_pre)
    trs = darkchem.utils.beamsearch(t_pred, k=k)
    trs = trs.reshape(-1,100)
    v2s = [darkchem.utils.vec2struct(trs[i]) for i in range(len(trs))]
    std = [standardize_smi(v2s[i]) for i in range(len(v2s))]
    return std


def pred_multiple_prod(testdf, model, path_vec, k=1):
    """Predicts the products of specific chemical reactions.
    Input is reactant SMILES strings.
    The default predicted consequence is one, you can change the value of k to get more probable forecasted results.

    Args:

    """
    a = []
    b = []
    c = []
    for i in range(len(testdf)):
        smi = testdf['Reactants'][i]
        std = transform_r2p_str(smi,model,path_vec,k)
        c.append(std)
        [a.append(std[i]) for i in range(len(std))]
    for j in range(len(std)):
        col = 'Product'
        b.append(col)
    out = pd.DataFrame(data = c, columns = b)
    out.insert(0,'Reactants',testdf['Reactants'].values,)
    df = struc2mol(pd.DataFrame(data = a,columns = ['smiles']))
    display(PandasTools.FrameToGridImage(df,column='mol', legendsCol='smiles',molsPerRow=5))
    return out


def pred_single_prod(smi,model,path_vec,k=1):
    """A function used to predict the product of a specific chemical reactions with the input of reactant smiles string.
    The default predicted consequence is one, you can change the value of k to get more probable forecasted results.
    """
    c = []
    b = []
    std = transform_r2p_str(smi, model, path_vec, k)
    c.append(std)
    for j in range(len(std)):
        col = 'Product'
        b.append(col)
    out = pd.DataFrame(data = c, columns = b)
    out.insert(0,'Reactant',smi)
    df = struc2mol(pd.DataFrame(data = std,columns = ['smiles']))
    display(PandasTools.FrameToGridImage(df,column='mol', legendsCol='smiles',molsPerRow=5))
    return out


def output_multiple_prod(testdf, model, path_vec, k=15):
    """A function used to output the product of many specific chemical reactions with the input of reactant smiles strings.
    The default value for k is 15.
    """
    a = []
    b = []
    c = []
    for i in range(len(testdf)):
        smi = testdf['Reactants'][i]
        a.append('Reactant')
        c.append(smi)
        std = transform_r2p_str(smi,model,path_vec,k)
        for j in range(len(std)):
            if std[j] == smi.upper():
                prd = std[j]
                break
            elif smi.replace('#','') == std[j]:
                prd = std[j]
                break      
            else:
                prd = std[14]
        a.append('Product')
        c.append(prd)
        b.append(prd)
    out = pd.DataFrame(data = b, columns = ['Products'])
    out.insert(0,'Reactants',testdf['Reactants'].values,)
    df = struc2mol(pd.DataFrame(data = c,columns = ['smiles']))
    df.insert(3,'legend',a)
    display(PandasTools.FrameToGridImage(df,column='mol', legendsCol='legend',molsPerRow=2))
    return out


def output_single_prod(smi,model,path_vec,k=15):  
    """A function used to predict the product of a specific chemical reactions with the input of reactant smiles string.
     When using beamsearch, the value of k is 15.
    """
    a = ['Reactant','Product']
    b = []
    c = [smi]
    std = transform_r2p_str(smi,model,path_vec,k)
    for j in range(len(std)):
        if std[j] == smi.upper(): # still need some more work, not applied to all reactions
            prd = std[j]
            break
        elif smi.replace('#','') == std[j]:
            prd = std[j]
            break      
        else:
            prd = std[14]
    b.append(prd)
    c.append(prd)
    out = pd.DataFrame(data = b, columns = ['Product'])
    out.insert(0,'Reactant',smi)
    df = struc2mol(pd.DataFrame(data = c,columns = ['smiles']))
    df.insert(3,'legend',a)
    display(PandasTools.FrameToGridImage(df,column='mol', legendsCol='legend',molsPerRow=5))
    return out
