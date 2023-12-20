#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(
    name='EMinsight',
    version='21',
    packages=find_packages(),
    description='A system for parsing TFS EPU data structures',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Kyle Morris',
    author_email='kyle.morris@diamond.ac.uk',
    url='https://github.com/kylelmorris/EMinsight',
    install_requires=[
        'glom',
        'tqdm',
        'pandas',
        'numpy',
        'scikit-learn',
        'matplotlib',
        'seaborn',
        'scipy',
        'xmltodict',
        'rich',
        'starparser',
        'mrcfile',
        'gemmi',
        'starfile',
        'pyem',
        'fpdf',
        'Pillow',
        'PyQt5',
        # Add any additional dependencies here
    ],
    python_requires='>=3.6', # Specify your Python version requirement
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        # More classifiers: https://pypi.org/classifiers/
    ],
)