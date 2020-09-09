#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_11_API_VERSION
#include <numpy/arrayobject.h>
#include "tbabsSrc/tbvabs.h"
#include "util.h"

static char tbabscore_docstring[] = 
    "This module wraps the tbabs XSPEC X-ray absorption package.";
static char tbnew_docstring[] = 
    "Compute X-ray absorption for arbitrary abundances following "
    "Wilms, Allen, and McCray, 2000, ApJ 542, 914-924";
static char vernabs_docstring[] = 
    "Compute cross section per Verner et al.";
static char dgami_docstring[] = 
    "Lower Incomplete Gamma Function reimplementation.";
static char phfit2_docstring[] = 
    "Cross sections directly from Verner et al 1996";
static char fgabnd_docstring[] = 
    "Abundances from Wilms+2000.";
static char deplf_docstring[] = 
    "Depletion factors from Wilms+2000.";

static PyObject *error_out(PyObject *m);
static PyObject *tbabscore_tbnew(PyObject *self, PyObject *args, 
                                 PyObject *kwargs);
static PyObject *tbabscore_vernabs(PyObject *self, PyObject *args, 
                                   PyObject *kwargs);
static PyObject *tbabscore_dgami(PyObject *self, PyObject *args, 
                                 PyObject *kwargs);
static PyObject *tbabscore_phfit2(PyObject *self, PyObject *args, 
                                 PyObject *kwargs);
//static PyObject *tbabscore_fgabnd(PyObject *self, PyObject *args, 
//                                 PyObject *kwargs);
//static PyObject *tbabscore_deplf(PyObject *self, PyObject *args, 
//                                 PyObject *kwargs);

// Declaring interface to Verner Fortran function.
void phfit2_(int *z,int *ion,int *shell,float *energ, float *sigma);

struct module_state
{
    PyObject *error;
};
#define GETSTATE(m) ((struct module_state *) PyModule_GetState(m))

static PyMethodDef tbabscoreMethods[] = {
    {"tbnew", (PyCFunction)tbabscore_tbnew, METH_VARARGS|METH_KEYWORDS,
        tbnew_docstring},
    {"vernabs", (PyCFunction)tbabscore_vernabs, METH_VARARGS|METH_KEYWORDS,
        vernabs_docstring},
    {"dgami", (PyCFunction)tbabscore_dgami, METH_VARARGS|METH_KEYWORDS,
        dgami_docstring},
    {"phfit2", (PyCFunction)tbabscore_phfit2, METH_VARARGS|METH_KEYWORDS,
        phfit2_docstring},
//    {"fgabnd", (PyCFunction)tbabscore_fgabnd, METH_VARARGS|METH_KEYWORDS,
//        fgabnd_docstring},
//    {"deplf", (PyCFunction)tbabscore_deplf, METH_VARARGS|METH_KEYWORDS,
//        deplf_docstring},
    {"error_out", (PyCFunction)error_out, METH_NOARGS, NULL},
    {NULL, NULL, 0, NULL}};


static int tbabscore_traverse(PyObject *m, visitproc visit, void *arg)
{
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int tbabscore_clear(PyObject *m)
{
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef tbabscoreModule = {
    PyModuleDef_HEAD_INIT,
    "tbabscore", /* Module Name */
    tbabscore_docstring,
    sizeof(struct module_state),
    tbabscoreMethods,
    NULL,
    tbabscore_traverse,
    tbabscore_clear,
    NULL
};
#define INITERROR return NULL

PyMODINIT_FUNC PyInit_jet(void)
{
    PyObject *module = PyModule_Create(&tbabscoreModule);
    if(module == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(module);
    st->error = PyErr_NewException("tbabscore.Error", NULL, NULL);
    if(st->error == NULL)
    {
        Py_DECREF(module);
        INITERROR;
    }

    //Load numpy stuff!
    import_array();

    return module;
}

static PyObject *error_out(PyObject *m)
{
    struct module_state *st = GETSTATE(m);
    PyErr_SetString(st->error, "something bad happened");
    return NULL;
}

static PyObject *tbabscore_tbnew(PyObject *self, PyObject *args, 
                                    PyObject *kwargs)
{
    PyObject *e_obj = NULL;
    PyObject *param_obj = NULL;

    static char *kwlist[] = {"Ebins", "param", NULL};

    //Parse Arguments
    if(!PyArg_ParseTupleAndKeywords(args, kwargs, "OO", kwlist,
                                    &e_obj, &param_obj))
    {
        PyErr_SetString(PyExc_RuntimeError, "Could not parse arguments.");
        return NULL;
    }

    //Grab NUMPY arrays
    PyArrayObject *e_arr;
    PyArrayObject *param_arr;

    e_arr = (PyArrayObject *) PyArray_FROM_OTF(e_obj, NPY_DOUBLE,
                                                NPY_ARRAY_IN_ARRAY);
    param_arr = (PyArrayObject *) PyArray_FROM_OTF(param_obj, NPY_DOUBLE,
                                                   NPY_ARRAY_IN_ARRAY);

    if(e_arr == NULL || param_arr == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, "Could not read input arrays.");
        Py_XDECREF(e_arr);
        Py_XDECREF(param_arr);
        return NULL;
    }

    // Check input array dimensions

    int e_ndim = (int) PyArray_NDIM(e_arr);
    int param_ndim = (int) PyArray_NDIM(param_arr);

    if(e_ndim != 1)
    {
        PyErr_SetString(PyExc_RuntimeError, "Ebins must be 1-D array");
        Py_DECREF(e_arr);
        Py_DECREF(param_arr);
        return NULL;
    }
    if(param_ndim != 1)
    {
        PyErr_SetString(PyExc_RuntimeError, "params must be 1-D array");
        Py_DECREF(e_arr);
        Py_DECREF(param_arr);
        return NULL;
    }


    // Check input array sizes
    int ne = (int)PyArray_DIM(e_arr, 0) - 1;
    int nparams = (int)PyArray_DIM(param_arr, 0);

    if(ne < 1)
    {
        PyErr_SetString(PyExc_RuntimeError, "Ebins must have length > 1");
        Py_DECREF(e_arr);
        Py_DECREF(param_arr);
        return NULL;
    }

    if(nparam != 42)
    {
        PyErr_SetString(PyExc_RuntimeError, "param must have length 42");
        Py_DECREF(e_arr);
        Py_DECREF(param_arr);
        return NULL;
    }

    double *e = (double *)PyArray_DATA(e_arr);
    double *param = (double *)PyArray_DATA(param_arr);

    //Allocate output array

    npy_intp dims[1] = {ne};
    PyObject *photar_obj = PyArray_SimpleNew(1, dims, NPY_DOUBLE);

    if(Fnu_obj == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, "Could not make photar array.");
        Py_DECREF(e_arr);
        Py_DECREF(param_arr);
        return NULL;
    }
    double *photar = PyArray_DATA((PyArrayObject *) photar_obj);

    // Calculate the absorption!
    tbnew(e, ne, param, 0, photarr, NULL, NULL);

    // Clean up!
    Py_DECREF(e_arr);
    Py_DECREF(param_arr);

    //Build output
    PyObject *ret = Py_BuildValue("N", photar_obj);
    
    return ret;
}

static PyObject *tbabscore_vernabs(PyObject *self, PyObject *args, 
                                    PyObject *kwargs)
{
    double e;
    int z;

    static char *kwlist[] = {"energy", "Z", NULL};

    //Parse Arguments
    if(!PyArg_ParseTupleAndKeywords(args, kwargs, "di", kwlist,
                                    &e, &z))
    {
        PyErr_SetString(PyExc_RuntimeError, "Could not parse arguments.");
        return NULL;
    }

    // Calculate the absorption!
    double sigma = vernabs(e, Z);

    //Build output
    PyObject *ret = Py_BuildValue("d", sigma);
    
    return ret;
}

static PyObject *tbabscore_phfit2(PyObject *self, PyObject *args, 
                                  PyObject *kwargs)
{
    int z, i, s;
    double e;

    static char *kwlist[] = {"Z", "I", "shell", "E", NULL};

    //Parse Arguments
    if(!PyArg_ParseTupleAndKeywords(args, kwargs, "iiid", kwlist,
                                    &z, &i, &s, &e))
    {
        PyErr_SetString(PyExc_RuntimeError, "Could not parse arguments.");
        return NULL;
    }

    // Calculate the absorption!
    float fsigma;
    float fe = (float)e;
    phfit2_(&z, &i, &s, &fe, &fsigma);
    double sigma = fsigma;

    //Build output
    PyObject *ret = Py_BuildValue("d", sigma);
    
    return ret;
}

static PyObject *tbabscore_dgami(PyObject *self, PyObject *args, 
                                 PyObject *kwargs)
{
    double s, x;

    static char *kwlist[] = {"s", "x", NULL};

    //Parse Arguments
    if(!PyArg_ParseTupleAndKeywords(args, kwargs, "dd", kwlist,
                                    &s, &x))
    {
        PyErr_SetString(PyExc_RuntimeError, "Could not parse arguments.");
        return NULL;
    }

    // Calculate the absorption!
    double y = gami(s, x);

    //Build output
    PyObject *ret = Py_BuildValue("d", y);
    
    return ret;
}
