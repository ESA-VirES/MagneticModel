#-------------------------------------------------------------------------------
#
#  data files
#
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
# The above copyright notice and this permission notice shall be included in all
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

from os.path import dirname, join

_DIRNAME = dirname(__file__)

# magnetic models
WMM_2015 = join(_DIRNAME, 'WMM2015v2.COF')
EMM_2010_STATIC = join(_DIRNAME, 'EMM-720_V3p0_static.cof')
EMM_2010_SECVAR = join(_DIRNAME, 'EMM-720_V3p0_secvar.cof')
CHAOS7_CORE_X9 = join(_DIRNAME, 'CHAOS-7.9_core.shc')
CHAOS7_CORE_X9_PREDICTION = join(_DIRNAME, 'CHAOS-7.9_core_prediction.shc')
CHAOS7_CORE_LATEST = CHAOS7_CORE_X9
CHAOS7_CORE_PREDICTION_LATEST = CHAOS7_CORE_X9_PREDICTION
CHAOS7_STATIC = join(_DIRNAME, 'CHAOS-7_static.shc')
CHAOS_STATIC_LATEST = CHAOS7_STATIC
CHAOS_CORE_LATEST = CHAOS7_CORE_LATEST
CHAOS_CORE_PREDICTION_LATEST = CHAOS7_CORE_PREDICTION_LATEST
IGRF11 = join(_DIRNAME, 'igrf11coeffs.txt')
IGRF12 = join(_DIRNAME, 'IGRF12.shc')
IGRF13 = join(_DIRNAME, 'IGRF13.shc')
IGRF14 = join(_DIRNAME, 'IGRF14.shc')
IGRF_LAST = IGRF14
IGRF_LAST_VERSION = "14"
SIFM = join(_DIRNAME, 'SIFM.shc')
LCS1 = join(_DIRNAME, 'LCS-1.shc')
MF7 = join(_DIRNAME, 'MF7.shc')

# magnetic models used by the apex point calculation
APEX_2015 = join(_DIRNAME, 'apexsh_1995-2015.txt')
APEX_2020 = join(_DIRNAME, 'apexsh_1980-2020.txt')
APEX_2025 = join(_DIRNAME, 'apexsh_1980-2025.txt')
APEX_2030 = join(_DIRNAME, 'apexsh_1980-2030.txt')
APEX_LAST = APEX_2030
