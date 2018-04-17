/*-----------------------------------------------------------------------------
 *
 * Time conversion utilities.
 *
 * Author: Martin Paces <martin.paces@eox.at>
 *
 *-----------------------------------------------------------------------------
 * Copyright (C) 2018 EOX IT Services GmbH
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
//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#define VERSION "0.1.0"

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

#include <string.h>
#include <stdlib.h>

/*---------------------------------------------------------------------------*/

#include "pytimeconv_mjd2000_to_decimal_year.h"
#include "pytimeconv_mjd2000_to_year_fraction.h"
#include "pytimeconv_decimal_year_to_mjd2000.h"

/*---------------------------------------------------------------------------*/
/* module's doc string */

#define DOC_PYSUNPOS \
"Time conversion utilities."

/*---------------------------------------------------------------------------*/
/*define module's methods */
static PyMethodDef pytimeconv_methods[] =
{
    {"mjd2000_to_decimal_year", (PyCFunction)pytimeconv_mjd2000_to_decimal_year, METH_VARARGS|METH_KEYWORDS, DOC_MJD2000_TO_DECIMAL_YEAR},
    {"mjd2000_to_year_fraction", (PyCFunction)pytimeconv_mjd2000_to_year_fraction, METH_VARARGS|METH_KEYWORDS, DOC_MJD2000_TO_YEAR_FRACTION},
    {"decimal_year_to_mjd2000", (PyCFunction)pytimeconv_decimal_year_to_mjd2000, METH_VARARGS|METH_KEYWORDS, DOC_DECIMAL_YEAR_TO_MJD2000},
    {NULL, NULL, 0, NULL} /* Sentinel - DO NOT REMOVE! */
} ;

/*---------------------------------------------------------------------------*/

/* module initialization  */
PyMODINIT_FUNC init_pytimeconv(void)
{
    PyObject *dict, *module;

    /* define module */
    module = Py_InitModule3("_pytimeconv", pytimeconv_methods, DOC_PYSUNPOS);
    if (NULL == module) return ;

    /* initialize Numpy arrays */
    import_array();

    dict = PyModule_GetDict(module);
    if (NULL == dict) return;

    /* metadata */
    PyDict_SetItemString(dict, "__author__", PyString_FromString("Martin Paces (martin.paces@eox.at)"));
    PyDict_SetItemString(dict, "__copyright__", PyString_FromString("Copyright (C) 2018 EOX IT Services GmbH"));
    PyDict_SetItemString(dict, "__licence__", PyString_FromString("EOX licence (MIT style)"));
    PyDict_SetItemString(dict, "__version__", PyString_FromString(VERSION));
}

/*---------------------------------------------------------------------------*/
