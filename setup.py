#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


<<<<<<< HEAD
with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

pkgs = find_packages(exclude=('example', 'docs', 'database'))
=======
# if os.path.exists('README.rst'):
#     long_description = open('README.rst').read()
# else:
#     long_description = '''An extension for DarkChem which predicts the products of chemical reactions and evaluates their encoding in DarkChem's latent space to return resultant properties.'''
>>>>>>> f47fa81cc30a21fcac53358e04246c99c99c6016

long_description = '''An extension for DarkChem which predicts the products of chemical reactions and evaluates their encoding in DarkChem's latent space to return resultant properties.'''
with open('LICENSE') as f:
    license = f.read()
    
setup(
    name='DarKnight',
    version='0.1',
    author='Christine Chang, Chih-Wei Hsu, Liang Xu',
<<<<<<< HEAD
    author_email='changch@uw.edu, cwh32@uw.edu, xuliang1@uw.edu',
    license=license,
=======
    author_email='changch@uw.edu, cwh32@uw.edu, xuliang1@uw.edu'
    license= license,
>>>>>>> f47fa81cc30a21fcac53358e04246c99c99c6016
    url='https://github.com/UWDIRECT-2019/DarKnight',
    packages=pkgs,
    description='Package to predict chemical reaction products and properties',
    long_description=readme,
    keywords='chemistry vae variational autoencoder latent space',
    # There is a problem for the rdkit package, we can't find its source from our terminal. So I deleted it temporarily.
    install_requires=[
        'darkchem', 'keras', 'tensorflow', 'openbabel', 'scikit-learn', 'mordred', 'numpy', 'scipy', 'matplotlib', 'seaborn', 'pandas'
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