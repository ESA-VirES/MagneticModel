#-------------------------------------------------------------------------------
#
#  Distutils Setup Script
#
# Project: Earth magnetic field in Python.
# Author: Martin Paces <martin.paces@eox.at>
#
#-------------------------------------------------------------------------------
# Copyright (C) 2014 EOX IT Services GmbH
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all
# copies of this Software or works derived from this Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#-------------------------------------------------------------------------------

import sys
from os.path import join
from distutils.core import setup
from distutils.extension import Extension
import eoxmagmod


COMMON_INCLUDE_DIRS = [
    './eoxmagmod',
    './eoxmagmod/include',
    join(sys.prefix, 'include'),
]

try:
    import numpy
    COMMON_INCLUDE_DIRS.append(numpy.get_include())
except ImportError:
    pass

setup(
    name="eoxmagmod",
    description="Earth magnetic field utilities.",
    author="Martin Paces",
    author_email="martin.paces@eox.at",
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: Implementation :: CPython',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Utilities',
    ],
    install_requires=[
        'numpy>=1.13.0',
        'spacepy',
    ],
    packages=[
        'eoxmagmod',
        'eoxmagmod.data',
        'eoxmagmod.tests',
        'eoxmagmod.tests.data',
        'eoxmagmod.magnetic_model',
        'eoxmagmod.magnetic_model.tests',
        'eoxmagmod.magnetic_model.tests.data',
    ],
    license='EOX licence (MIT style)',
    version=eoxmagmod.__version__,
    package_data={
        'eoxmagmod': [
            'data/*',
            'tests/data/*.tsv',
            'magnetic_model/tests/data/*.txt',
            'magnetic_model/tests/data/*.cdf',
        ],
    },
    ext_modules=[
        Extension(
            'eoxmagmod._pymm',
            sources=[
                'eoxmagmod/pymm.c',
            ],
            libraries=[],
            library_dirs=[],
            include_dirs=COMMON_INCLUDE_DIRS,
        ),
        Extension(
            'eoxmagmod._pyqd',
            sources=[
                'eoxmagmod/pyqd.c',
            ],
            libraries=['qdipole'],
            library_dirs=[],
            include_dirs=COMMON_INCLUDE_DIRS,
        ),
        Extension(
            'eoxmagmod._pysunpos',
            sources=[
                'eoxmagmod/pysunpos.c',
            ],
            libraries=[],
            library_dirs=[],
            include_dirs=COMMON_INCLUDE_DIRS,
        ),
        Extension(
            'eoxmagmod._pytimeconv',
            sources=[
                'eoxmagmod/pytimeconv.c',
            ],
            libraries=[],
            library_dirs=[],
            include_dirs=COMMON_INCLUDE_DIRS,
        ),
    ]
)
