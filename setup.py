#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from setuptools import setup, find_packages


if os.path.exists('README.rst'):
    long_description = open('README.rst').read()
else:
    long_description = '''An extension for DarkChem which predicts the products of chemical reactions and evaluates their encoding in DarkChem's latent space to return resultant properties.'''

with open('LICENSE') as f:
    license = f.read()
    
setup(
    name='DarKnight',
    version='0.1',
    author='Christine Chang, Chih-Wei Hsu, Liang Xu',
    author_email='changch@uw.edu, cwh32@uw.edu, xuliang1@uw.edu'
    license='MIT',
    url='https://github.com/UWDIRECT-2019/DarKnight',
    packages=find_packages(),
    description='Package to predict chemical reaction products and properties',
    long_description=long_description,
    keywords='chemistry vae variational autoencoder latent space',
    install_requires=[
        'darkchem', 'rdkit', 'keras', 'tensorflow', 'openbabel', 'scikit-learn', 'mordred', 'numpy', 'scipy', 'matplotlib', 'seaborn', 'pandas'
    ],
    classifiers=[
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Chemistry',
    ],
)
