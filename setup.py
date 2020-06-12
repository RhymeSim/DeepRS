"""PyRhyme setup file"""

import os
from setuptools import setup


def read(fname):
    """Reading a given file with its relative path"""
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name='ml_riemann_solver',
    version='2020.6',
    author='Saeed Sarpas',
    author_email='saeed.sarpas@phys.ethz.ch',
    description='Machine Learning Powered Riemann Solver',
    long_description=read('README.md'),
    long_description_content_type="text/markdown",
    url='git@gitlab.com:rhyme-org/ml_riemann_solver.git',
    keywords=['Rhyme', 'Machine Learning',  'Riemann Solver'],
    license='GPLv3',
    packages=['ml_riemann_solver'],
    install_requires=[
        'numpy',
    ],
    scripts=[],
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    zip_safe=False,
)
