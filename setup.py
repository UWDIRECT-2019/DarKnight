#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

pkgs = find_packages(exclude=('example', 'docs', 'database'))

setup(
    name='DarKnight',
    version='0.1',
    author='Christine Chang, Chih-Wei Hsu, Liang Xu',
    author_email='changch@uw.edu, cwh32@uw.edu, xuliang1@uw.edu',
    license=license,
    url='https://github.com/UWDIRECT-2019/DarKnight',
    packages=pkgs,
    description='Package to predict chemical reaction products and properties',
    long_description=readme,
    keywords='chemistry vae variational autoencoder latent space',
    install_requires=[
        'darkchem', 'keras', 'tensorflow','rdkit', 'openbabel', 'scikit-learn', 'mordred', 'numpy', 'scipy', 'matplotlib', 'seaborn', 'pandas'
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
