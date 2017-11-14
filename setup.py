# -*- coding: utf-8 -*-
import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "antiana",
    version = "0.0.1",
    author = "Jiaming Meng",
    author_email = "jiamingm1990@outlook.com",
    description = ("antibody analysis tools on protein sequencing"),
    license = "Apache License",
    keywords = "antibody analysis",
    url = "https://github.com/jeremy1990/AntibodyAnalysis",
    packages = ['antiana', 'docs'],
    long_description = read('README.rst'),
    classifiers = [
        "Development Status :: 1 - Planning",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: Apache Software License"
    ]
)
