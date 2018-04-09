#-------------------------------------------------------------------------------
#
#  Earth Magnetic Model 2010
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

import re
from gzip import GzipFile as _GzipFile
from numpy import zeros
from .wmm import MagneticModelSimple
from .data import EMM_2010_STATIC, EMM_2010_SECVAR


# NOTE: GzipFile with statement support added in Python 2.7.
if hasattr(_GzipFile, '__exit__'):
    GzipFile = _GzipFile
else:
    class GzipFile(_GzipFile):
        def __enter__(self):
            return self
        def __exit__(self, type_, value, traceback):
            self.close()


def _read_prm_emm2010(fname):
    """ read single coefficient file """

    def _read(fid):
        prm = {'sources': [fname]}
        #PASS 1
        # parse header
        headers = []
        line = next(fid)
        name = line.split()[0]
        if name[-1] == ',':
            name = name[:-1]
        prm['name'] = name
        headers.append(line.rstrip())
        lidx = 0
        for _lidx, line in enumerate(fid, 1):
            lidx = _lidx
            #print "%d\t%s"%(lidx, line.rstrip())
            headers.append(line.rstrip())
            wlist = line.split()
            if wlist and wlist[0] == 'Epoch:':
                epoch = float(wlist[1])
            elif wlist and wlist[0] == 'Degree:':
                degree = int(wlist[1])
            if re.match(r'^\s*-+\s*$', line):
                break

        # parse coefficients
        prm['headers'] = ["\n".join(headers)]
        prm['degree'] = degree
        prm['epoch'] = epoch

        degree_check = 0
        nterm = ((degree+2)*(degree+1))/2
        coef_g = zeros(nterm)
        coef_h = zeros(nterm)

        loff = lidx
        for lidx, line in enumerate(fid):
            #print loff+lidx, line.rstrip()
            try:
                n, m, g, h = line.split()
                n, m = int(n), int(m)
                g, h, = float(g), float(h)
                if m > n:
                    raise ValueError
            except ValueError:
                raise ValueError(
                    "Invalid line #%d: %s" % (loff + lidx, line.rstrip())
                )
            else:
                degree_check = max(degree_check, n)
                if degree_check > degree:
                    raise ValueError(
                        "Invalid order value! The annotated order value is too "
                        "low!"
                    )

                idx = m + ((n+1)*n)/2
                coef_g[idx] = g
                coef_h[idx] = h
                #lcoef.append((n, m, g, h, dg, dh))

        if degree_check > degree:
            raise ValueError(
                "Invalid order value! The annotated order value is too high!"
            )

        prm['coef_g'] = coef_g
        prm['coef_h'] = coef_h

        return prm

    try:
        with GzipFile(fname, 'r') as fid:
            return _read(fid)
    except IOError:
        with file(fname, 'r') as fid:
            return _read(fid)


def read_model_emm2010(fname_static=EMM_2010_STATIC,
                       fname_secvar=EMM_2010_SECVAR):
    """ Read model parameters from the coefficient files in the WMM2010 format.
    """

    prm_static = _read_prm_emm2010(fname_static)
    prm_secvar = _read_prm_emm2010(fname_secvar)

    if prm_static['name'] != prm_secvar['name']:
        raise ValueError("Model name mismatch in the coefficients' files!")
    if prm_static['epoch'] != prm_secvar['epoch']:
        raise ValueError("Model epoch mismatch in the coefficients' files!")

    prm = {}
    prm['name'] = prm_static['name']
    prm['version'] = prm_static['name']
    prm['epoch'] = prm_static['epoch']
    prm['valid_until'] = prm_static['epoch'] + 5.0
    prm['sources'] = prm_static['sources'] + prm_secvar['sources']
    prm['headers'] = prm_static['headers'] + prm_secvar['headers']
    prm['degree_static'] = prm_static['degree']
    prm['degree_secvar'] = prm_secvar['degree']
    prm['coef_static_g'] = prm_static['coef_g']
    prm['coef_static_h'] = prm_static['coef_h']
    prm['coef_secvar_g'] = prm_secvar['coef_g']
    prm['coef_secvar_h'] = prm_secvar['coef_h']

    return MagneticModelSimple(prm)
