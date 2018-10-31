/*-----------------------------------------------------------------------------
 *
 * Magnetic Quasi Dipole Coordinates - C python bindings
 *
 * Author: Martin Paces <martin.paces@eox.at>
 *
 *-----------------------------------------------------------------------------
 * Copyright (C) 2015 EOX IT Services GmbH
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies of this Software or works derived from this Software.
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

/* module specific exceptions */
#include "pyqd_exc.h"

/* quasi dipole coordinates evaluation - old API with wrong MLT */
#include "pyqd_eval_apex.h"

/* Quasi-Dipole coordinates evaluation */
#include "pyqd_eval_qdlatlon.h"

/* Magnetic Local Time coordinates evaluation */
#include "pyqd_eval_mlt.h"

/* sub-solar point coordinates evaluation */
#include "pyqd_eval_subsol.h"
/*---------------------------------------------------------------------------*/
/* module's doc string */

#define DOC_PYQD \
"library evaluation magnetic quasi-dipole coordinates and local magnetic time."

/*---------------------------------------------------------------------------*/
/*define module's methods */
static PyMethodDef pyqd_methods[] =
{
    {"eval_qdlatlon", (PyCFunction)eval_qdlatlon, METH_VARARGS|METH_KEYWORDS, DOC_EVAL_QDLATLON},
    {"eval_mlt", (PyCFunction)eval_mlt, METH_VARARGS|METH_KEYWORDS, DOC_EVAL_MLT},
    {"eval_subsol", (PyCFunction)eval_subsol, METH_VARARGS|METH_KEYWORDS, DOC_EVAL_SUBSOL},
    {"eval_apex", (PyCFunction)eval_apex, METH_VARARGS|METH_KEYWORDS, DOC_EVAL_APEX},
    {NULL, NULL, 0, NULL} /* Sentinel - DO NOT REMOVE! */
} ;

/*---------------------------------------------------------------------------*/

/* module initialization  */
PyMODINIT_FUNC init_pyqd(void)
{
    PyObject *dict, *module;

    /* define module */
    if (NULL == (module = Py_InitModule3("_pyqd", pyqd_methods, DOC_PYQD)))
        return;

    /* initialize numpy arrays */
    import_array();

    if (NULL == (dict = PyModule_GetDict(module)))
        return;

    /* add RT2Error */
    PyExc_QDError = PyErr_NewException("_pyqd.QDError", NULL, NULL);

    PyDict_SetItemString(dict, "QDError", PyExc_QDError);

    /* constants */
    #define SET_INT_ITEM(d, s, i) \
    {PyObject *tmp = PyInt_FromLong(i);  PyDict_SetItemString(d,s,tmp); Py_DECREF(tmp);}

    /* metadata */
    PyDict_SetItemString(dict, "__author__", PyString_FromString("Martin Paces (martin.paces@eox.at)"));
    PyDict_SetItemString(dict, "__copyright__", PyString_FromString("Copyright (C) 2015 EOX IT Services GmbH"));
    PyDict_SetItemString(dict, "__licence__", PyString_FromString("EOX licence (MIT style)"));
    PyDict_SetItemString(dict, "__version__", PyString_FromString(VERSION));
}

/*---------------------------------------------------------------------------*/
