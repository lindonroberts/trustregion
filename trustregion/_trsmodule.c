/* File: _trsmodule.c
 * This file is auto-generated with f2py (version:1.26.4).
 * f2py is a Fortran to Python Interface Generator (FPIG), Second Edition,
 * written by Pearu Peterson <pearu@cens.ioc.ee>.
 * Generation date: Mon Apr 15 04:53:11 2024
 * Do not edit this file directly unless you know what you are doing!!!
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef PY_SSIZE_T_CLEAN
#define PY_SSIZE_T_CLEAN
#endif /* PY_SSIZE_T_CLEAN */

/* Unconditionally included */
#include <Python.h>
#include <numpy/npy_os.h>

/*********************** See f2py2e/cfuncs.py: includes ***********************/
#include "fortranobject.h"
#include <math.h>

/**************** See f2py2e/rules.py: mod_rules['modulebody'] ****************/
static PyObject *_trs_error;
static PyObject *_trs_module;

/*********************** See f2py2e/cfuncs.py: typedefs ***********************/
/*need_typedefs*/

/****************** See f2py2e/cfuncs.py: typedefs_generated ******************/
/*need_typedefs_generated*/

/********************** See f2py2e/cfuncs.py: cppmacros **********************/

#if defined(PREPEND_FORTRAN)
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F
#else
#define F_FUNC(f,F) _##f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F##_
#else
#define F_FUNC(f,F) _##f##_
#endif
#endif
#else
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F
#else
#define F_FUNC(f,F) f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F##_
#else
#define F_FUNC(f,F) f##_
#endif
#endif
#endif
#if defined(UNDERSCORE_G77)
#define F_FUNC_US(f,F) F_FUNC(f##_,F##_)
#else
#define F_FUNC_US(f,F) F_FUNC(f,F)
#endif


/* See fortranobject.h for definitions. The macros here are provided for BC. */
#define rank f2py_rank
#define shape f2py_shape
#define fshape f2py_shape
#define len f2py_len
#define flen f2py_flen
#define slen f2py_slen
#define size f2py_size


#ifdef DEBUGCFUNCS
#define CFUNCSMESS(mess) fprintf(stderr,"debug-capi:"mess);
#define CFUNCSMESSPY(mess,obj) CFUNCSMESS(mess) \
    PyObject_Print((PyObject *)obj,stderr,Py_PRINT_RAW);\
    fprintf(stderr,"\n");
#else
#define CFUNCSMESS(mess)
#define CFUNCSMESSPY(mess,obj)
#endif


#ifndef max
#define max(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef min
#define min(a,b) ((a < b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) ((a < b) ? (a) : (b))
#endif


/************************ See f2py2e/cfuncs.py: cfuncs ************************/

static int
double_from_pyobj(double* v, PyObject *obj, const char *errmess)
{
    PyObject* tmp = NULL;
    if (PyFloat_Check(obj)) {
        *v = PyFloat_AsDouble(obj);
        return !(*v == -1.0 && PyErr_Occurred());
    }

    tmp = PyNumber_Float(obj);
    if (tmp) {
        *v = PyFloat_AsDouble(tmp);
        Py_DECREF(tmp);
        return !(*v == -1.0 && PyErr_Occurred());
    }

    if (PyComplex_Check(obj)) {
        PyErr_Clear();
        tmp = PyObject_GetAttrString(obj,"real");
    }
    else if (PyBytes_Check(obj) || PyUnicode_Check(obj)) {
        /*pass*/;
    }
    else if (PySequence_Check(obj)) {
        PyErr_Clear();
        tmp = PySequence_GetItem(obj, 0);
    }

    if (tmp) {
        if (double_from_pyobj(v,tmp,errmess)) {Py_DECREF(tmp); return 1;}
        Py_DECREF(tmp);
    }
    {
        PyObject* err = PyErr_Occurred();
        if (err==NULL) err = _trs_error;
        PyErr_SetString(err,errmess);
    }
    return 0;
}


/********************* See f2py2e/cfuncs.py: userincludes *********************/
/*need_userincludes*/

/********************* See f2py2e/capi_rules.py: usercode *********************/


/* See f2py2e/rules.py */
extern void F_FUNC(trsapp,TRSAPP)(int*,double*,double*,double*,double*,double*,double*);
extern void F_FUNC(trsbox,TRSBOX)(int*,double*,double*,double*,double*,double*,double*,double*,double*,double*);
extern void F_FUNC(trslin,TRSLIN)(int*,double*,double*,double*,double*,double*,double*);
/*eof externroutines*/

/******************** See f2py2e/capi_rules.py: usercode1 ********************/


/******************* See f2py2e/cb_rules.py: buildcallback *******************/
/*need_callbacks*/

/*********************** See f2py2e/rules.py: buildapi ***********************/

/*********************************** trsapp ***********************************/
static char doc_f2py_rout__trs_trsapp[] = "\
STEP,CRVMIN = trsapp(XOPT,GQ,HQ,DELTA)\n\nWrapper for ``trsapp``.\
\n\nParameters\n----------\n"
"XOPT : input rank-1 array('d') with bounds (N)\n"
"GQ : input rank-1 array('d') with bounds (N)\n"
"HQ : input rank-2 array('d') with bounds (N,N)\n"
"DELTA : input float\n"
"\nReturns\n-------\n"
"STEP : rank-1 array('d') with bounds (N)\n"
"CRVMIN : float";
/* extern void F_FUNC(trsapp,TRSAPP)(int*,double*,double*,double*,double*,double*,double*); */
static PyObject *f2py_rout__trs_trsapp(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(int*,double*,double*,double*,double*,double*,double*)) {
    PyObject * volatile capi_buildvalue = NULL;
    volatile int f2py_success = 1;
/*decl*/

    int N = 0;
    double *XOPT = NULL;
    npy_intp XOPT_Dims[1] = {-1};
    const int XOPT_Rank = 1;
    PyArrayObject *capi_XOPT_as_array = NULL;
    int capi_XOPT_intent = 0;
    PyObject *XOPT_capi = Py_None;
    double *GQ = NULL;
    npy_intp GQ_Dims[1] = {-1};
    const int GQ_Rank = 1;
    PyArrayObject *capi_GQ_as_array = NULL;
    int capi_GQ_intent = 0;
    PyObject *GQ_capi = Py_None;
    double *HQ = NULL;
    npy_intp HQ_Dims[2] = {-1, -1};
    const int HQ_Rank = 2;
    PyArrayObject *capi_HQ_as_array = NULL;
    int capi_HQ_intent = 0;
    PyObject *HQ_capi = Py_None;
    double DELTA = 0;
    PyObject *DELTA_capi = Py_None;
    double *STEP = NULL;
    npy_intp STEP_Dims[1] = {-1};
    const int STEP_Rank = 1;
    PyArrayObject *capi_STEP_as_array = NULL;
    int capi_STEP_intent = 0;
    double CRVMIN = 0;
    static char *capi_kwlist[] = {"XOPT","GQ","HQ","DELTA",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
    if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
        "OOOO|:_trs.trsapp",\
        capi_kwlist,&XOPT_capi,&GQ_capi,&HQ_capi,&DELTA_capi))
        return NULL;
/*frompyobj*/
    /* Processing variable DELTA */
        f2py_success = double_from_pyobj(&DELTA,DELTA_capi,"_trs.trsapp() 4th argument (DELTA) can't be converted to double");
    if (f2py_success) {
    /* Processing variable XOPT */
    ;
    capi_XOPT_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "_trs._trs.trsapp: failed to create array from the 1st argument `XOPT`";
    capi_XOPT_as_array = ndarray_from_pyobj(  NPY_DOUBLE,1,XOPT_Dims,XOPT_Rank,  capi_XOPT_intent,XOPT_capi,capi_errmess);
    if (capi_XOPT_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = _trs_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        XOPT = (double *)(PyArray_DATA(capi_XOPT_as_array));

    /* Processing variable CRVMIN */
    /* Processing variable N */
    N = shape(XOPT,0);
    /* Processing variable GQ */
    GQ_Dims[0]=N;
    capi_GQ_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "_trs._trs.trsapp: failed to create array from the 2nd argument `GQ`";
    capi_GQ_as_array = ndarray_from_pyobj(  NPY_DOUBLE,1,GQ_Dims,GQ_Rank,  capi_GQ_intent,GQ_capi,capi_errmess);
    if (capi_GQ_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = _trs_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        GQ = (double *)(PyArray_DATA(capi_GQ_as_array));

    /* Processing variable HQ */
    HQ_Dims[0]=N,HQ_Dims[1]=N;
    capi_HQ_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "_trs._trs.trsapp: failed to create array from the 3rd argument `HQ`";
    capi_HQ_as_array = ndarray_from_pyobj(  NPY_DOUBLE,1,HQ_Dims,HQ_Rank,  capi_HQ_intent,HQ_capi,capi_errmess);
    if (capi_HQ_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = _trs_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        HQ = (double *)(PyArray_DATA(capi_HQ_as_array));

    /* Processing variable STEP */
    STEP_Dims[0]=N;
    capi_STEP_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
    const char * capi_errmess = "_trs._trs.trsapp: failed to create array from the hidden `STEP`";
    capi_STEP_as_array = ndarray_from_pyobj(  NPY_DOUBLE,1,STEP_Dims,STEP_Rank,  capi_STEP_intent,Py_None,capi_errmess);
    if (capi_STEP_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = _trs_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        STEP = (double *)(PyArray_DATA(capi_STEP_as_array));

/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
                (*f2py_func)(&N,XOPT,GQ,HQ,&DELTA,STEP,&CRVMIN);
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
        if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
        CFUNCSMESS("Building return value.\n");
        capi_buildvalue = Py_BuildValue("Nd",capi_STEP_as_array,CRVMIN);
/*closepyobjfrom*/
/*end of closepyobjfrom*/
        } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
    }  /* if (capi_STEP_as_array == NULL) ... else of STEP */
    /* End of cleaning variable STEP */
    if((PyObject *)capi_HQ_as_array!=HQ_capi) {
        Py_XDECREF(capi_HQ_as_array); }
    }  /* if (capi_HQ_as_array == NULL) ... else of HQ */
    /* End of cleaning variable HQ */
    if((PyObject *)capi_GQ_as_array!=GQ_capi) {
        Py_XDECREF(capi_GQ_as_array); }
    }  /* if (capi_GQ_as_array == NULL) ... else of GQ */
    /* End of cleaning variable GQ */
    /* End of cleaning variable N */
    /* End of cleaning variable CRVMIN */
    if((PyObject *)capi_XOPT_as_array!=XOPT_capi) {
        Py_XDECREF(capi_XOPT_as_array); }
    }  /* if (capi_XOPT_as_array == NULL) ... else of XOPT */
    /* End of cleaning variable XOPT */
    } /*if (f2py_success) of DELTA*/
    /* End of cleaning variable DELTA */
/*end of cleanupfrompyobj*/
    if (capi_buildvalue == NULL) {
/*routdebugfailure*/
    } else {
/*routdebugleave*/
    }
    CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
    return capi_buildvalue;
}
/******************************* end of trsapp *******************************/

/*********************************** trsbox ***********************************/
static char doc_f2py_rout__trs_trsbox[] = "\
D,GNEW,CRVMIN = trsbox(XOPT,GOPT,HQ,SL,SU,DELTA)\n\nWrapper for ``trsbox``.\
\n\nParameters\n----------\n"
"XOPT : input rank-1 array('d') with bounds (N)\n"
"GOPT : input rank-1 array('d') with bounds (N)\n"
"HQ : input rank-2 array('d') with bounds (N,N)\n"
"SL : input rank-1 array('d') with bounds (N)\n"
"SU : input rank-1 array('d') with bounds (N)\n"
"DELTA : input float\n"
"\nReturns\n-------\n"
"D : rank-1 array('d') with bounds (N)\n"
"GNEW : rank-1 array('d') with bounds (N)\n"
"CRVMIN : float";
/* extern void F_FUNC(trsbox,TRSBOX)(int*,double*,double*,double*,double*,double*,double*,double*,double*,double*); */
static PyObject *f2py_rout__trs_trsbox(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(int*,double*,double*,double*,double*,double*,double*,double*,double*,double*)) {
    PyObject * volatile capi_buildvalue = NULL;
    volatile int f2py_success = 1;
/*decl*/

    int N = 0;
    double *XOPT = NULL;
    npy_intp XOPT_Dims[1] = {-1};
    const int XOPT_Rank = 1;
    PyArrayObject *capi_XOPT_as_array = NULL;
    int capi_XOPT_intent = 0;
    PyObject *XOPT_capi = Py_None;
    double *GOPT = NULL;
    npy_intp GOPT_Dims[1] = {-1};
    const int GOPT_Rank = 1;
    PyArrayObject *capi_GOPT_as_array = NULL;
    int capi_GOPT_intent = 0;
    PyObject *GOPT_capi = Py_None;
    double *HQ = NULL;
    npy_intp HQ_Dims[2] = {-1, -1};
    const int HQ_Rank = 2;
    PyArrayObject *capi_HQ_as_array = NULL;
    int capi_HQ_intent = 0;
    PyObject *HQ_capi = Py_None;
    double *SL = NULL;
    npy_intp SL_Dims[1] = {-1};
    const int SL_Rank = 1;
    PyArrayObject *capi_SL_as_array = NULL;
    int capi_SL_intent = 0;
    PyObject *SL_capi = Py_None;
    double *SU = NULL;
    npy_intp SU_Dims[1] = {-1};
    const int SU_Rank = 1;
    PyArrayObject *capi_SU_as_array = NULL;
    int capi_SU_intent = 0;
    PyObject *SU_capi = Py_None;
    double DELTA = 0;
    PyObject *DELTA_capi = Py_None;
    double *D = NULL;
    npy_intp D_Dims[1] = {-1};
    const int D_Rank = 1;
    PyArrayObject *capi_D_as_array = NULL;
    int capi_D_intent = 0;
    double *GNEW = NULL;
    npy_intp GNEW_Dims[1] = {-1};
    const int GNEW_Rank = 1;
    PyArrayObject *capi_GNEW_as_array = NULL;
    int capi_GNEW_intent = 0;
    double CRVMIN = 0;
    static char *capi_kwlist[] = {"XOPT","GOPT","HQ","SL","SU","DELTA",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
    if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
        "OOOOOO|:_trs.trsbox",\
        capi_kwlist,&XOPT_capi,&GOPT_capi,&HQ_capi,&SL_capi,&SU_capi,&DELTA_capi))
        return NULL;
/*frompyobj*/
    /* Processing variable DELTA */
        f2py_success = double_from_pyobj(&DELTA,DELTA_capi,"_trs.trsbox() 6th argument (DELTA) can't be converted to double");
    if (f2py_success) {
    /* Processing variable XOPT */
    ;
    capi_XOPT_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "_trs._trs.trsbox: failed to create array from the 1st argument `XOPT`";
    capi_XOPT_as_array = ndarray_from_pyobj(  NPY_DOUBLE,1,XOPT_Dims,XOPT_Rank,  capi_XOPT_intent,XOPT_capi,capi_errmess);
    if (capi_XOPT_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = _trs_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        XOPT = (double *)(PyArray_DATA(capi_XOPT_as_array));

    /* Processing variable CRVMIN */
    /* Processing variable N */
    N = shape(XOPT,0);
    /* Processing variable GOPT */
    GOPT_Dims[0]=N;
    capi_GOPT_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "_trs._trs.trsbox: failed to create array from the 2nd argument `GOPT`";
    capi_GOPT_as_array = ndarray_from_pyobj(  NPY_DOUBLE,1,GOPT_Dims,GOPT_Rank,  capi_GOPT_intent,GOPT_capi,capi_errmess);
    if (capi_GOPT_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = _trs_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        GOPT = (double *)(PyArray_DATA(capi_GOPT_as_array));

    /* Processing variable SL */
    SL_Dims[0]=N;
    capi_SL_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "_trs._trs.trsbox: failed to create array from the 4th argument `SL`";
    capi_SL_as_array = ndarray_from_pyobj(  NPY_DOUBLE,1,SL_Dims,SL_Rank,  capi_SL_intent,SL_capi,capi_errmess);
    if (capi_SL_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = _trs_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        SL = (double *)(PyArray_DATA(capi_SL_as_array));

    /* Processing variable SU */
    SU_Dims[0]=N;
    capi_SU_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "_trs._trs.trsbox: failed to create array from the 5th argument `SU`";
    capi_SU_as_array = ndarray_from_pyobj(  NPY_DOUBLE,1,SU_Dims,SU_Rank,  capi_SU_intent,SU_capi,capi_errmess);
    if (capi_SU_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = _trs_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        SU = (double *)(PyArray_DATA(capi_SU_as_array));

    /* Processing variable HQ */
    HQ_Dims[0]=N,HQ_Dims[1]=N;
    capi_HQ_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "_trs._trs.trsbox: failed to create array from the 3rd argument `HQ`";
    capi_HQ_as_array = ndarray_from_pyobj(  NPY_DOUBLE,1,HQ_Dims,HQ_Rank,  capi_HQ_intent,HQ_capi,capi_errmess);
    if (capi_HQ_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = _trs_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        HQ = (double *)(PyArray_DATA(capi_HQ_as_array));

    /* Processing variable D */
    D_Dims[0]=N;
    capi_D_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
    const char * capi_errmess = "_trs._trs.trsbox: failed to create array from the hidden `D`";
    capi_D_as_array = ndarray_from_pyobj(  NPY_DOUBLE,1,D_Dims,D_Rank,  capi_D_intent,Py_None,capi_errmess);
    if (capi_D_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = _trs_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        D = (double *)(PyArray_DATA(capi_D_as_array));

    /* Processing variable GNEW */
    GNEW_Dims[0]=N;
    capi_GNEW_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
    const char * capi_errmess = "_trs._trs.trsbox: failed to create array from the hidden `GNEW`";
    capi_GNEW_as_array = ndarray_from_pyobj(  NPY_DOUBLE,1,GNEW_Dims,GNEW_Rank,  capi_GNEW_intent,Py_None,capi_errmess);
    if (capi_GNEW_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = _trs_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        GNEW = (double *)(PyArray_DATA(capi_GNEW_as_array));

/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
                (*f2py_func)(&N,XOPT,GOPT,HQ,SL,SU,&DELTA,D,GNEW,&CRVMIN);
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
        if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
        CFUNCSMESS("Building return value.\n");
        capi_buildvalue = Py_BuildValue("NNd",capi_D_as_array,capi_GNEW_as_array,CRVMIN);
/*closepyobjfrom*/
/*end of closepyobjfrom*/
        } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
    }  /* if (capi_GNEW_as_array == NULL) ... else of GNEW */
    /* End of cleaning variable GNEW */
    }  /* if (capi_D_as_array == NULL) ... else of D */
    /* End of cleaning variable D */
    if((PyObject *)capi_HQ_as_array!=HQ_capi) {
        Py_XDECREF(capi_HQ_as_array); }
    }  /* if (capi_HQ_as_array == NULL) ... else of HQ */
    /* End of cleaning variable HQ */
    if((PyObject *)capi_SU_as_array!=SU_capi) {
        Py_XDECREF(capi_SU_as_array); }
    }  /* if (capi_SU_as_array == NULL) ... else of SU */
    /* End of cleaning variable SU */
    if((PyObject *)capi_SL_as_array!=SL_capi) {
        Py_XDECREF(capi_SL_as_array); }
    }  /* if (capi_SL_as_array == NULL) ... else of SL */
    /* End of cleaning variable SL */
    if((PyObject *)capi_GOPT_as_array!=GOPT_capi) {
        Py_XDECREF(capi_GOPT_as_array); }
    }  /* if (capi_GOPT_as_array == NULL) ... else of GOPT */
    /* End of cleaning variable GOPT */
    /* End of cleaning variable N */
    /* End of cleaning variable CRVMIN */
    if((PyObject *)capi_XOPT_as_array!=XOPT_capi) {
        Py_XDECREF(capi_XOPT_as_array); }
    }  /* if (capi_XOPT_as_array == NULL) ... else of XOPT */
    /* End of cleaning variable XOPT */
    } /*if (f2py_success) of DELTA*/
    /* End of cleaning variable DELTA */
/*end of cleanupfrompyobj*/
    if (capi_buildvalue == NULL) {
/*routdebugfailure*/
    } else {
/*routdebugleave*/
    }
    CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
    return capi_buildvalue;
}
/******************************* end of trsbox *******************************/

/*********************************** trslin ***********************************/
static char doc_f2py_rout__trs_trslin[] = "\
D,CRVMIN = trslin(GOPT,SL,SU,DELTA)\n\nWrapper for ``trslin``.\
\n\nParameters\n----------\n"
"GOPT : input rank-1 array('d') with bounds (N)\n"
"SL : input rank-1 array('d') with bounds (N)\n"
"SU : input rank-1 array('d') with bounds (N)\n"
"DELTA : input float\n"
"\nReturns\n-------\n"
"D : rank-1 array('d') with bounds (N)\n"
"CRVMIN : float";
/* extern void F_FUNC(trslin,TRSLIN)(int*,double*,double*,double*,double*,double*,double*); */
static PyObject *f2py_rout__trs_trslin(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(int*,double*,double*,double*,double*,double*,double*)) {
    PyObject * volatile capi_buildvalue = NULL;
    volatile int f2py_success = 1;
/*decl*/

    int N = 0;
    double *GOPT = NULL;
    npy_intp GOPT_Dims[1] = {-1};
    const int GOPT_Rank = 1;
    PyArrayObject *capi_GOPT_as_array = NULL;
    int capi_GOPT_intent = 0;
    PyObject *GOPT_capi = Py_None;
    double *SL = NULL;
    npy_intp SL_Dims[1] = {-1};
    const int SL_Rank = 1;
    PyArrayObject *capi_SL_as_array = NULL;
    int capi_SL_intent = 0;
    PyObject *SL_capi = Py_None;
    double *SU = NULL;
    npy_intp SU_Dims[1] = {-1};
    const int SU_Rank = 1;
    PyArrayObject *capi_SU_as_array = NULL;
    int capi_SU_intent = 0;
    PyObject *SU_capi = Py_None;
    double DELTA = 0;
    PyObject *DELTA_capi = Py_None;
    double *D = NULL;
    npy_intp D_Dims[1] = {-1};
    const int D_Rank = 1;
    PyArrayObject *capi_D_as_array = NULL;
    int capi_D_intent = 0;
    double CRVMIN = 0;
    static char *capi_kwlist[] = {"GOPT","SL","SU","DELTA",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
    if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
        "OOOO|:_trs.trslin",\
        capi_kwlist,&GOPT_capi,&SL_capi,&SU_capi,&DELTA_capi))
        return NULL;
/*frompyobj*/
    /* Processing variable DELTA */
        f2py_success = double_from_pyobj(&DELTA,DELTA_capi,"_trs.trslin() 4th argument (DELTA) can't be converted to double");
    if (f2py_success) {
    /* Processing variable GOPT */
    ;
    capi_GOPT_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "_trs._trs.trslin: failed to create array from the 1st argument `GOPT`";
    capi_GOPT_as_array = ndarray_from_pyobj(  NPY_DOUBLE,1,GOPT_Dims,GOPT_Rank,  capi_GOPT_intent,GOPT_capi,capi_errmess);
    if (capi_GOPT_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = _trs_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        GOPT = (double *)(PyArray_DATA(capi_GOPT_as_array));

    /* Processing variable CRVMIN */
    /* Processing variable N */
    N = shape(GOPT,0);
    /* Processing variable SL */
    SL_Dims[0]=N;
    capi_SL_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "_trs._trs.trslin: failed to create array from the 2nd argument `SL`";
    capi_SL_as_array = ndarray_from_pyobj(  NPY_DOUBLE,1,SL_Dims,SL_Rank,  capi_SL_intent,SL_capi,capi_errmess);
    if (capi_SL_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = _trs_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        SL = (double *)(PyArray_DATA(capi_SL_as_array));

    /* Processing variable SU */
    SU_Dims[0]=N;
    capi_SU_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "_trs._trs.trslin: failed to create array from the 3rd argument `SU`";
    capi_SU_as_array = ndarray_from_pyobj(  NPY_DOUBLE,1,SU_Dims,SU_Rank,  capi_SU_intent,SU_capi,capi_errmess);
    if (capi_SU_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = _trs_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        SU = (double *)(PyArray_DATA(capi_SU_as_array));

    /* Processing variable D */
    D_Dims[0]=N;
    capi_D_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
    const char * capi_errmess = "_trs._trs.trslin: failed to create array from the hidden `D`";
    capi_D_as_array = ndarray_from_pyobj(  NPY_DOUBLE,1,D_Dims,D_Rank,  capi_D_intent,Py_None,capi_errmess);
    if (capi_D_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = _trs_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        D = (double *)(PyArray_DATA(capi_D_as_array));

/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
                (*f2py_func)(&N,GOPT,SL,SU,&DELTA,D,&CRVMIN);
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
        if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
        CFUNCSMESS("Building return value.\n");
        capi_buildvalue = Py_BuildValue("Nd",capi_D_as_array,CRVMIN);
/*closepyobjfrom*/
/*end of closepyobjfrom*/
        } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
    }  /* if (capi_D_as_array == NULL) ... else of D */
    /* End of cleaning variable D */
    if((PyObject *)capi_SU_as_array!=SU_capi) {
        Py_XDECREF(capi_SU_as_array); }
    }  /* if (capi_SU_as_array == NULL) ... else of SU */
    /* End of cleaning variable SU */
    if((PyObject *)capi_SL_as_array!=SL_capi) {
        Py_XDECREF(capi_SL_as_array); }
    }  /* if (capi_SL_as_array == NULL) ... else of SL */
    /* End of cleaning variable SL */
    /* End of cleaning variable N */
    /* End of cleaning variable CRVMIN */
    if((PyObject *)capi_GOPT_as_array!=GOPT_capi) {
        Py_XDECREF(capi_GOPT_as_array); }
    }  /* if (capi_GOPT_as_array == NULL) ... else of GOPT */
    /* End of cleaning variable GOPT */
    } /*if (f2py_success) of DELTA*/
    /* End of cleaning variable DELTA */
/*end of cleanupfrompyobj*/
    if (capi_buildvalue == NULL) {
/*routdebugfailure*/
    } else {
/*routdebugleave*/
    }
    CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
    return capi_buildvalue;
}
/******************************* end of trslin *******************************/
/*eof body*/

/******************* See f2py2e/f90mod_rules.py: buildhooks *******************/
/*need_f90modhooks*/

/************** See f2py2e/rules.py: module_rules['modulebody'] **************/

/******************* See f2py2e/common_rules.py: buildhooks *******************/

/*need_commonhooks*/

/**************************** See f2py2e/rules.py ****************************/

static FortranDataDef f2py_routine_defs[] = {
    {"trsapp",-1,{{-1}},0,0,(char *)  F_FUNC(trsapp,TRSAPP),  (f2py_init_func)f2py_rout__trs_trsapp,doc_f2py_rout__trs_trsapp},
    {"trsbox",-1,{{-1}},0,0,(char *)  F_FUNC(trsbox,TRSBOX),  (f2py_init_func)f2py_rout__trs_trsbox,doc_f2py_rout__trs_trsbox},
    {"trslin",-1,{{-1}},0,0,(char *)  F_FUNC(trslin,TRSLIN),  (f2py_init_func)f2py_rout__trs_trslin,doc_f2py_rout__trs_trslin},

/*eof routine_defs*/
    {NULL}
};

static PyMethodDef f2py_module_methods[] = {

    {NULL,NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_trs",
    NULL,
    -1,
    f2py_module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit__trs(void) {
    int i;
    PyObject *m,*d, *s, *tmp;
    m = _trs_module = PyModule_Create(&moduledef);
    Py_SET_TYPE(&PyFortran_Type, &PyType_Type);
    import_array();
    if (PyErr_Occurred())
        {PyErr_SetString(PyExc_ImportError, "can't initialize module _trs (failed to import numpy)"); return m;}
    d = PyModule_GetDict(m);
    s = PyUnicode_FromString("1.26.4");
    PyDict_SetItemString(d, "__version__", s);
    Py_DECREF(s);
    s = PyUnicode_FromString(
        "This module '_trs' is auto-generated with f2py (version:1.26.4).\nFunctions:\n"
"    STEP,CRVMIN = trsapp(XOPT,GQ,HQ,DELTA)\n"
"    D,GNEW,CRVMIN = trsbox(XOPT,GOPT,HQ,SL,SU,DELTA)\n"
"    D,CRVMIN = trslin(GOPT,SL,SU,DELTA)\n"
".");
    PyDict_SetItemString(d, "__doc__", s);
    Py_DECREF(s);
    s = PyUnicode_FromString("1.26.4");
    PyDict_SetItemString(d, "__f2py_numpy_version__", s);
    Py_DECREF(s);
    _trs_error = PyErr_NewException ("_trs.error", NULL, NULL);
    /*
     * Store the error object inside the dict, so that it could get deallocated.
     * (in practice, this is a module, so it likely will not and cannot.)
     */
    PyDict_SetItemString(d, "__trs_error", _trs_error);
    Py_DECREF(_trs_error);
    for(i=0;f2py_routine_defs[i].name!=NULL;i++) {
        tmp = PyFortranObject_NewAsAttr(&f2py_routine_defs[i]);
        PyDict_SetItemString(d, f2py_routine_defs[i].name, tmp);
        Py_DECREF(tmp);
    }



/*eof initf2pywraphooks*/
/*eof initf90modhooks*/

/*eof initcommonhooks*/


#ifdef F2PY_REPORT_ATEXIT
    if (! PyErr_Occurred())
        on_exit(f2py_report_on_exit,(void*)"_trs");
#endif
    return m;
}
#ifdef __cplusplus
}
#endif