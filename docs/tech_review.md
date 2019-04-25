# Technology Review

DarkChem uses Python 3.6.  The required dependencies are:

- `openbabel`
- `keras`
- `tensorflow`
- `rdkit`
- `openbabel`
- `numpy`
- `scipy`
- `scikit-learn`
- `matplotlib` 
- `seaborn` 
- `pandas` 

## Packages Used in DarkChem

#### `multiprocessing`:Support spawning processes package 
##### https://docs.python.org/3.4/library/multiprocessing.html?highlight=process

#### `re`: `Regular Expression Operations`: Strings checking package 
##### https://docs.python.org/2/library/re.html

#### `glob`: Documenets managing and pathnames finding module 
##### https://docs.python.org/2/library/glob.html

#### `pandas`: Python Data Analysis Library
##### https://pandas.pydata.org/

#### `subprocess`: This package alloes you to spawn new processes 
##### https://docs.python.org/2/library/subprocess.html 

#### `os`: Miscellaneous operating system interfaces 
##### https://docs.python.org/3/library/os.html

#### `keras`: Neural-network and Deep Learning library 
##### https://keras.io/

#### `numpy`: Calculation package for arrays and matrixs
##### http://cs231n.github.io/python-numpy-tutorial/

#### `ast`: `Abstract Syntax Trees`: Pcakage helping programmers to parse and revise pyhton codes 
##### https://docs.python.org/3/library/ast.html

#### `rdkit (Chem)`: Transforming SMILES into Python coding 
##### https://www.rdkit.org/docs/GettingStartedInPython.html

#### `functools (partial)`: Module permiting the immobilization of certain parameters of a function, then the function is renewed. 
##### https://docs.python.org/3/library/functools.html

#### `argparse`: A Python built in package, which is utilized to analyze options and parameters of command lines 
##### https://docs.python.org/3/library/argparse.html

# Understanding DarkChem

## General

DarkChem uses a **variational autoencoder (VAE)** in order to learn a latent (numerical vector space) representation of chemicals, given their SMILES string representations.  The VAE includes an [encoder layer](https://www.quora.com/What-is-an-Encoder-Decoder-in-Deep-Learning), which encodes each molecule as a vector in this latent space.  The VAE also includes a decoder layer, to translate the latent space representation into desired physical properties (_m/z_ ratio and collision cross-section (CCS)).  The decoder can map these properties onto the raw format (the initial molecule).

[Transfer learning](https://towardsdatascience.com/a-comprehensive-hands-on-guide-to-transfer-learning-with-real-world-applications-in-deep-learning-212bf3b2f27a) was used to develop the model which eventually produced the highest *reconstruction accuracy* (what does this mean).  In the case of DarkChem, transfer learning proceeded through the following iterations upon training networks:

    N1b > N3b > N7b
    
The `N7b` trained models are what were eventually provided to us by PNNL for our implementation.

## Networks

### N1b

Properties of the N1b network:

* **Mode:** VAE, _m/z_
* **Data:** PubChem
* **Weights:** None
* **Accuracy**
  * Reconstruction
  * _m/z_ Error
  * CCS Error
  
### N3b

* **Mode:** VAE, _m/z_, collision x-section (CCS)
* **Data:** _in silico_
* **Weights:** inherited from N1b
* **Accuracy**
  * Reconstruction
  * _m/z_ Error
  * CCS Error
  
### N7b

* **Mode:** VAE, _m/z_, collision x-section (CCS)
* **Data:** experimental
* **Weights:** inherited from N3b, with weights frozen
* **Accuracy**
  * Reconstruction
  * _m/z_ Error
  * CCS Error