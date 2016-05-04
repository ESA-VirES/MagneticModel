/*-----------------------------------------------------------------------------
 *
 * World Magnetic Model - C python bindings
 *
 * Project: World Magnetic Model - python interface
 * Author: Martin Paces <martin.paces@eox.at>
 *
 *-----------------------------------------------------------------------------
 * Copyright (C) 2014 EOX IT Services GmbH
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

#define VERSION "0.2.0"

// needed to prevent dual definition
#ifdef _POSIX_C_SOURCE
#undef _POSIX_C_SOURCE
#endif

#define PY_SSIZE_T_CLEAN 1
#include <Python.h>
#include <numpy/arrayobject.h>

#if PY_MAJOR_VERSION != 2
#error "Non-supported Python major version!"
#endif
#if PY_MINOR_VERSION < 6
#error "Non-supported Python minor version!"
#endif

/* maximum allowed output array dimension */
#define MAX_OUT_ARRAY_NDIM 16

/*---------------------------------------------------------------------------*/

/* Module specific exceptions. */
#include "pywmm_exc.h"

/* Coordinate conversions. */
#include "pywmm_cconv.h"

/* spherical harmonic model evaluation */
#include "pywmm_sheval.h"

/* evaluation of the associative Legendre functions */
#include "pywmm_legendre.h"

/* evaluation of the relative radius powers */
#include "pywmm_relradpow.h"

/* evaluation of the series of longitude sine and cosine values */
#include "pywmm_lonsincos.h"

/* final spherical-harmonic gradient evaluation */
#include "pywmm_sphargrd.h"

/* final spherical-harmonic gradient potential */
#include "pywmm_spharpot.h"

/* vector rotations */
#include "pywmm_vrot_sph2geod.h"
#include "pywmm_vrot_sph2cart.h"
#include "pywmm_vrot_cart2sph.h"

/*---------------------------------------------------------------------------*/
/* module's doc string */

#define DOC_PYWMM \
"This module provides bindings to the World Magnetic Model library."

/*---------------------------------------------------------------------------*/
/*define module's methods */
static PyMethodDef pywmm_methods[] =
{
    {"vrot_sph2cart", (PyCFunction)vrot_sph2cart, METH_VARARGS|METH_KEYWORDS, DOC_VROT_SPH2CART},
    {"vrot_cart2sph", (PyCFunction)vrot_cart2sph, METH_VARARGS|METH_KEYWORDS, DOC_VROT_CART2SPH},
    {"vrot_sph2geod", (PyCFunction)vrot_sph2geod, METH_VARARGS|METH_KEYWORDS, DOC_VROT_SPH2GEOD},
    {"spharpot", (PyCFunction)spharpot, METH_VARARGS|METH_KEYWORDS, DOC_SPHARPOT},
    {"sphargrd", (PyCFunction)sphargrd, METH_VARARGS|METH_KEYWORDS, DOC_SPHARGRD},
    {"lonsincos", (PyCFunction)lonsincos, METH_VARARGS|METH_KEYWORDS, DOC_LONSINCOS},
    {"relradpow", (PyCFunction)relradpow, METH_VARARGS|METH_KEYWORDS, DOC_RELRADPOW},
    {"legendre", (PyCFunction)legendre, METH_VARARGS|METH_KEYWORDS, DOC_LEGENDRE},
    {"sheval", (PyCFunction)sheval, METH_VARARGS|METH_KEYWORDS, DOC_SHEVAL},
    {"convert", (PyCFunction)convert, METH_VARARGS|METH_KEYWORDS, DOC_CONVERT},
    {NULL, NULL, 0, NULL} /* Sentinel - DO NOT REMOVE! */
} ;

/*---------------------------------------------------------------------------*/

/* module initialization  */
PyMODINIT_FUNC init_pywmm(void)
{
    PyObject *dict, *module;

    /* define module */
    module = Py_InitModule3("_pywmm", pywmm_methods, DOC_PYWMM);
    if (NULL == module) return ;

    /* initialize numpy arrays */
    import_array();

    dict = PyModule_GetDict(module);
    if (NULL == dict) return;

    /* add RT2Error */
    PyExc_WMMError = PyErr_NewException("_pywmm.WMMError", NULL, NULL);

    PyDict_SetItemString(dict, "WMMError", PyExc_WMMError);

    /* constants */
    #define SET_INT_ITEM(d, s, i) \
    {PyObject *tmp = PyInt_FromLong(i);  PyDict_SetItemString(d,s,tmp); Py_DECREF(tmp);}
    SET_INT_ITEM(dict, "GEODETIC_ABOVE_WGS84", CT_GEODETIC_ABOVE_WGS84);
    SET_INT_ITEM(dict, "GEODETIC_ABOVE_EGM96", CT_GEODETIC_ABOVE_EGM96);
    SET_INT_ITEM(dict, "GEOCENTRIC_SPHERICAL", CT_GEOCENTRIC_SPHERICAL);
    SET_INT_ITEM(dict, "GEOCENTRIC_CARTESIAN", CT_GEOCENTRIC_CARTESIAN);
    SET_INT_ITEM(dict, "POTENTIAL", SM_POTENTIAL);
    SET_INT_ITEM(dict, "GRADIENT", SM_GRADIENT);
    SET_INT_ITEM(dict, "POTENTIAL_AND_GRADIENT", SM_POTENTIAL_AND_GRADIENT);

    /* metadata */
    PyDict_SetItemString(dict, "__author__", PyString_FromString("Martin Paces (martin.paces@eox.at)"));
    PyDict_SetItemString(dict, "__copyright__", PyString_FromString("Copyright (C) 2014 EOX IT Services GmbH"));
    PyDict_SetItemString(dict, "__licence__", PyString_FromString("EOX licence (MIT style)"));
    PyDict_SetItemString(dict, "__version__", PyString_FromString(VERSION));
}

/*---------------------------------------------------------------------------*/
