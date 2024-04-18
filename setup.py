#!/usr/bin/env python
# -*- coding: utf-8 -*-

# {# pkglts, pysetup.kwds
# format setup arguments

from os import walk
from os.path import abspath, normpath
from os.path import join as pj

from setuptools import setup, find_packages


short_descr = "ArchiCrop is a parametric plant model that span cereals and legume"
#readme = open('README.rst').read()
#history = open('HISTORY.rst').read().replace('.. :changelog:', '')
readme = short_descr
history = ''
# find version number in src/hydroshoot/version.py
version = {}
with open("src/archicrop/version.py") as fp:
    exec(fp.read(), version)

'''
data_files = []

nb = len(normpath(abspath("src/hydroshoot_data"))) + 1


def data_rel_pth(pth):
    """ Return path relative to pkg_data
    """
    abs_pth = normpath(abspath(pth))
    return abs_pth[nb:]


for root, dnames, fnames in walk("src/hydroshoot_data"):
    for name in fnames:
        data_files.append(data_rel_pth(pj(root, name)))
'''

setup_kwds = dict(
    name='openalea.archicrop',
    version=version["__version__"],
    description=short_descr,
    long_description=readme + '\n\n' + history,
    author="Oriane & co ",
    author_email="oriane braud at cirad francais",
    url='https://github.com/openalea/archicrop',
    license='cecill-c',
    zip_safe=False,

    packages=find_packages('src'),
    package_dir={'': 'src'},
    
 #   include_package_data=True,
#    package_data={'hydroshoot_data': data_files},
    install_requires=[
        ],
    tests_require=[
        "pytest",
        ],
    entry_points={},
    keywords='',
    test_suite='nose.collector',
)
# #}
# change setup_kwds below before the next pkglts tag

# do not change things below
# {# pkglts, pysetup.call
setup(**setup_kwds)
# #}
