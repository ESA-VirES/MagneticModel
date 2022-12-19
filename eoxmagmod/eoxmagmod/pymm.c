/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings
 *
 * Author: Martin Paces <martin.paces@eox.at>
 *
 *-----------------------------------------------------------------------------
 * Copyright (C) 2014-2022 EOX IT Services GmbH
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies of this Software or works derived from this Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *-----------------------------------------------------------------------------
*/

#include "common.h" /* common definitions - to be included before Python.h */
#include <Python.h>
#include <numpy/arrayobject.h>

/* module version */
#include "version.h"

/* common python utilities */
#include "py_common.h"

/* Coordinate conversions. */
#include "pymm_cconv.h"

/* spherical harmonic model evaluation */
#include "pymm_sheval.h"
#include "pymm_shevaltemp.h"
#include "pymm_sheval2dfs.h"

/* evaluation of the associative Legendre functions */
#include "pymm_legendre.h"

/* evaluation of the relative radius powers */
#include "pymm_relradpow.h"

/* evaluation of the series of longitude cosine and sine values */
#include "pymm_loncossin.h"

/* final spherical-harmonic gradient evaluation */
#include "pymm_sphargrd.h"

/* final spherical-harmonic gradient potential */
#include "pymm_spharpot.h"

/* vector rotations */
#include "pymm_vrot_sph2geod.h"
#include "pymm_vrot_sph2cart.h"
#include "pymm_vrot_cart2sph.h"

/* bisect interval search */
#include "pymm_bisect.h"

/* coefficients interpolation */
#include "pymm_interp.h"

/* 2D Fourier series coefficients evaluation */
#include "pymm_fourier2d.h"

/*---------------------------------------------------------------------------*/
/* module's doc string */

#define DOC_PYMM \
"This module provides bindings to the Geomagnetic Model library."

/*---------------------------------------------------------------------------*/
/*define module's methods */
static PyMethodDef pymm_methods[] =
{
    {"vrot_sph2cart", (PyCFunction)vrot_sph2cart, METH_VARARGS|METH_KEYWORDS, DOC_VROT_SPH2CART},
    {"vrot_cart2sph", (PyCFunction)vrot_cart2sph, METH_VARARGS|METH_KEYWORDS, DOC_VROT_CART2SPH},
    {"vrot_sph2geod", (PyCFunction)vrot_sph2geod, METH_VARARGS|METH_KEYWORDS, DOC_VROT_SPH2GEOD},
    {"spharpot", (PyCFunction)spharpot, METH_VARARGS|METH_KEYWORDS, DOC_SPHARPOT},
    {"sphargrd", (PyCFunction)sphargrd, METH_VARARGS|METH_KEYWORDS, DOC_SPHARGRD},
    {"loncossin", (PyCFunction)loncossin, METH_VARARGS|METH_KEYWORDS, DOC_LONCOSSIN},
    {"relradpow", (PyCFunction)relradpow, METH_VARARGS|METH_KEYWORDS, DOC_RELRADPOW},
    {"legendre", (PyCFunction)legendre, METH_VARARGS|METH_KEYWORDS, DOC_LEGENDRE},
    {"sheval", (PyCFunction)sheval, METH_VARARGS|METH_KEYWORDS, DOC_SHEVAL},
    {"shevaltemp", (PyCFunction)shevaltemp, METH_VARARGS|METH_KEYWORDS, DOC_SHEVALTEMP},
    {"sheval2dfs", (PyCFunction)sheval2dfs, METH_VARARGS|METH_KEYWORDS, DOC_SHEVAL2DFS},
    {"convert", (PyCFunction)convert, METH_VARARGS|METH_KEYWORDS, DOC_CONVERT},
    {"bisect", (PyCFunction)bisect, METH_VARARGS|METH_KEYWORDS, DOC_BISECT},
    {"interp", (PyCFunction)interp, METH_VARARGS|METH_KEYWORDS, DOC_INTERP},
    {"fourier2d", (PyCFunction)fourier2d, METH_VARARGS|METH_KEYWORDS, DOC_FOURIER2D},
    {NULL, NULL, 0, NULL} /* sentinel - DO NOT REMOVE! */
} ;

/*---------------------------------------------------------------------------*/
/* module initialization  */

static PyObject* init_module(void)
{
    PyObject *module = init_python_module("_pymm", DOC_PYMM, pymm_methods);
    if (NULL == module)
        goto exit;

    PyObject *dict = PyModule_GetDict(module);
    if (NULL == dict)
        goto exit;

    /* add module specific exception */
    //set_dict_item_str(dict, "MMError", PyErr_NewException("_pymm.MMError", NULL, NULL));

    /* integer constants */
    set_dict_item_str_long(dict, "GEODETIC_ABOVE_WGS84", CT_GEODETIC_ABOVE_WGS84);
    set_dict_item_str_long(dict, "GEOCENTRIC_SPHERICAL", CT_GEOCENTRIC_SPHERICAL);
    set_dict_item_str_long(dict, "GEOCENTRIC_CARTESIAN", CT_GEOCENTRIC_CARTESIAN);
    set_dict_item_str_long(dict, "POTENTIAL", SM_POTENTIAL);
    set_dict_item_str_long(dict, "GRADIENT", SM_GRADIENT);
    set_dict_item_str_long(dict, "POTENTIAL_AND_GRADIENT", SM_POTENTIAL_AND_GRADIENT);
    set_dict_item_str_double(dict, "EARTH_RADIUS", RADIUS);
    set_dict_item_str_long(dict, "BISECT_SIDE_LEFT", BISECT_SIDE_LEFT);
    set_dict_item_str_long(dict, "BISECT_SIDE_RIGHT", BISECT_SIDE_RIGHT);
    set_dict_item_str_long(dict, "INTERP_C0", INTERP_C0);
    set_dict_item_str_long(dict, "INTERP_C1", INTERP_C1);
    set_dict_item_str_long(dict, "INTERP_C1D1", INTERP_C1D1);

    /* module metadata */
    set_dict_item_str_str(dict, "__author__", "Martin Paces (martin.paces@eox.at)");
    set_dict_item_str_str(dict, "__copyright__", "Copyright (C) 2014-2022 EOX IT Services GmbH");
    set_dict_item_str_str(dict, "__licence__", "EOX licence (MIT style)");
    set_dict_item_str_str(dict, "__version__", VERSION);

  exit:
    return module;
}

/*---------------------------------------------------------------------------*/

PyObject* PyInit__pymm(void)
{
    import_array();
    return init_module();
}
