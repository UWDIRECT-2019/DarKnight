# DarKnight
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://travis-ci.org/UWDIRECT-2019/DarKnight.svg?branch=master)](https://travis-ci.org/UWDIRECT-2019/DarKnight)
[![Coverage Status](https://coveralls.io/repos/github/UWDIRECT-2019/DarKnight/badge.svg?branch=master)](https://coveralls.io/github/UWDIRECT-2019/DarKnight?branch=master)

## Overview
The goal of our project is to explore the potential use of a PNNL-developed software package, DarkChem, to relate chemical transformations to observed physical properties. DarkChem uses a variational autoencoder (VAE) to encode molecular structure and properties in latent (numerical vector) space; the reverse transformation (decoding) allows for the retrieval of structures that match desired chemical properties. In our capstone project, we will explore how two molecules in this latent space, related through chemical reactions, are related. Understanding and identifying this relationship enables the application of DarkChem towards the prediction of unknown products in a chemical reaction, as well as their molecular properties.

<div align=center> <img src="https://github.com/UWDIRECT-2019/DarKnight/raw/master/figures/logo.jpg" width="400"> </div>

## Use Cases
#### 1. Predict unknown product(s) in a given reaction, as well as the properties of the product(s).
Input(s): Molecular structure of reactant(s) (SMILES), type or class of reaction (e.g. oxidation, reduction, etc.)
Output(s): Molecular structure of product(s) (SMILES), properties of product(s) (e.g. m/z, CCS)
Components: Input structure, convert structure to latent space representation, perform reaction vector transformation on structure in latent space, use DarkChem to transform latent space representation into structure and property values
#### 2. Predict the structure of a small molecule given desired molecular properties.
Input(s): Desired properties (e.g. m/z, CCS)
Output(s): Molecular structure of product(s) matching the specified properties (SMILES, or visual representation)
Components: Input properties, convert properties into latent space representation using DarkChem, decode latent space representation to retrieve structures, convert latent space structure representation to SMILES or other visualization
#### 3. Predict properties for a new or unknown molecule given its structure.
Input(s): Molecular structure of the product (as a SMILES string)
Output(s): Molecular properties predicted for the structure (e.g. m/z, CCS)
Components: Input structure, convert structure to latent space representation, use DarkChem to relate latent space representation to properties, output properties
#### 4. Components: Input structure, convert structure to latent space representation, use DarkChem to relate latent space representation to properties, output properties
Input(s): Molecular structures for reactants (as SMILES strings)
Output(s): Molecular properties predicted for the structure (e.g. m/z, CCS)
Components: Input structure(s) and/or properties, convert structure(s) to latent space representations, use DarkChem to search latent space for potential transformations matching inputs, output possible product structures

## Package Requirements

* darkchem
* rdkit 
* keras 
* tensorflow
* openbabel
* scikit-learn
* mordred
* numpy
* scipy 
* matplotlib
* seaborn
* pandas

## Installation
------------
Use [``conda``](https://www.anaconda.com/download/) to create a new virtual environment with required dependencies from the `environment.yml`:
```bash
conda env create -f environment.yml
```

Activate the virtual environment:
```
conda activate darknight
```

Install DarkChem using [``pip``](https://pypi.org/project/pip/):
```bash
# direct
pip install git+https://github.com/pnnl/darkchem.git
```

## DarKnight Instruction
<div align=center> <img src="https://github.com/UWDIRECT-2019/DarKnight/raw/master/figures/GUI_gif.gif" width="400"> </div>
As the GIF showed above, what we need to do is input our reactant smile strings and select the exactly type of reaction, then click the Predict button, finally we will get the result in 15 seconds
