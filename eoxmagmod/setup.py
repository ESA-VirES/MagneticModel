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

from distutils.core import setup
from distutils.extension import Extension

setup(
    name="eoxmagmod",
    description="Earth magnetic field utilities.",
    author="Martin Paces",
    author_email="martin.paces@eox.at",
    packages=['eoxmagmod'],
    license='EOX licence (MIT style)',
    version='0.4.0',
    package_data={'eoxmagmod': ['data/*']},
    ext_modules=[
        Extension(
            'eoxmagmod._pywmm',
            sources=[
                'eoxmagmod/pywmm.c',
            ],
            libraries=['geomag'],
            library_dirs=[],
            include_dirs=[
                './eoxmagmod',
                './eoxmagmod/include',
            ],
        ),
        Extension(
            'eoxmagmod._pyqd',
            sources=[
                'eoxmagmod/pyqd.c',
            ],
            libraries=['qdipole'],
            library_dirs=[],
            include_dirs=[
                './eoxmagmod',
                './eoxmagmod/include',
            ],
        ),
        Extension(
            'eoxmagmod._pysunpos',
            sources=[
                'eoxmagmod/pysunpos.c',
            ],
            libraries=[],
            library_dirs=[],
            include_dirs=[
                './eoxmagmod',
                './eoxmagmod/include',
            ],
        ),
    ]
)
