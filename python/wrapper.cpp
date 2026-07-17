#define NPY_NO_DEPRECATED_API NPY_2_0_API_VERSION
#include <Python.h>
#include <memory>
#include "../src/driver/driver.h"
#include <numpy/arrayobject.h>  // 包含 NumPy 头文件

// 类型定义
typedef struct {
    PyObject_HEAD
    twinkle::driver_t* obj;  // 核心对象指针
} PyTwinkle;

// tp_new 实现（仅分配内存）
static PyObject* PyTwinkle_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    // 调用基类分配内存
    PyTwinkle *self = (PyTwinkle*)type->tp_alloc(type, 0);
    if (!self) return NULL;
    // 初始化为空指针，防止野指针
    self->obj = nullptr;
    return (PyObject*)self;
}

static void PyTwinkle_dealloc(PyTwinkle *self) {
    if (self->obj) {
        self->obj->free();   // 删除 source_base 内容
        delete self->obj;    // 触发 C++ 析构函数
        self->obj = nullptr; // 防止野指针
    }
    // printf("deleted!\n");
    Py_TYPE(self)->tp_free((PyObject*)self);
}

// init()方法包装
static int PyTwinkle_init(PyTwinkle* self, PyObject* args, PyObject* kwargs) {
    // 参数默认值
    int n_srcs = -1;
    int device_num = 0;
    int n_stream = 1;
    double RelTol = 1e-4;
    int astrom = 0;
    int pt_only = 0;

    // 定义参数选项
    static char* kwlist1[] = {(char*)"n_srcs", (char*)"device_num", (char*)"n_stream", (char*)"RelTol", NULL};
    static char* kwlist2[] = {(char*)"n_srcs", (char*)"device_num", (char*)"n_stream", NULL};
    static char* kwlist3[] = {(char*)"n_srcs", (char*)"device_num", (char*)"n_stream", (char*)"RelTol", (char*)"astrom", NULL};
    static char* kwlist4[] = {(char*)"n_srcs", (char*)"device_num", (char*)"n_stream", (char*)"RelTol", (char*)"astrom", (char*)"pt_only", NULL};
    
    // 解析参数字典
    if  (PyArg_ParseTupleAndKeywords(args, kwargs, "i|iidii", kwlist4, &n_srcs, &device_num, &n_stream, &RelTol, &astrom, &pt_only))
    {    }
    else
    {
        PyErr_Clear();
        if  (PyArg_ParseTupleAndKeywords(args, kwargs, "iii|di", kwlist3, &n_srcs, &device_num, &n_stream, &RelTol, &astrom))
        {    }
        else
        {
            PyErr_Clear();
            if (PyArg_ParseTupleAndKeywords(args, kwargs, "iii|d", kwlist1, &n_srcs, &device_num, &n_stream, &RelTol))
            {    }
            else
            {
                PyErr_Clear();
                if (PyArg_ParseTupleAndKeywords(args, kwargs, "iii", kwlist2, &n_srcs, &device_num, &n_stream))
                {  
                    // printf("In new version, RelTol is set in init(), default 1e-4. The RelTol in set_params() won't be used\n");
                    printf("RelTol should be set in init(), default 1e-4.\n");
                }
                else
                {
                    return -1;
                }
            }        
        }
    }



    try {
        self->obj = new twinkle::driver_t();
        self->obj->init(n_srcs, device_num, n_stream, RelTol, (bool)astrom, (bool)pt_only);
    } 
    catch (const std::bad_alloc& e) {
        PyErr_SetString(PyExc_MemoryError, "Memory allocation failed");
        return -1; // 必须返回-1
    }
    catch (const std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return -1;
    }
    catch (...) {
        PyErr_SetString(PyExc_RuntimeError, "Unknown exception");
        return -1;
    }

    return 0;
}


static PyObject* PyTwinkle_set_params(PyTwinkle *self, PyObject *args) {
    PyObject *py_ss = nullptr;
    double qq = 0.0, rho = 0.0, RelTol = -1.0;
    PyObject *py_xs = nullptr, *py_ys = nullptr;
    PyArrayObject *arr_xs = nullptr, *arr_ys = nullptr, *arr_ss = nullptr;
    npy_intp len_x = 0, len_y = 0; // 提前声明变量

    // 解析参数
    if (PyArg_ParseTuple(args, "OdddOO", &py_ss, &qq, &rho, &RelTol, &py_xs, &py_ys)) {
        // 第一种情况：提供了RelTol
    }
    else {
        PyErr_Clear();
        if (!PyArg_ParseTuple(args, "OddOO", &py_ss, &qq, &rho, &py_xs, &py_ys)) {
            PyErr_SetString(PyExc_TypeError, 
                "Arguments: (ss: float OR np.ndarray, qq, rho: float, "
                "xs: np.ndarray, ys: np.ndarray) OR "
                "(ss: float OR np.ndarray, qq, rho: float, RelTol: float, "
                "xs: np.ndarray, ys: np.ndarray)");
            return NULL;
        }
    }

    // 转换并检查 xs 数组
    arr_xs = (PyArrayObject*)PyArray_FromAny(py_xs, 
        PyArray_DescrFromType(NPY_DOUBLE), 1, 1,
        NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED, NULL);
    if (!arr_xs) {
        PyErr_SetString(PyExc_TypeError, "xs must be a 1D contiguous float64 array");
        goto cleanup;
    }

    // 转换并检查 ys 数组
    arr_ys = (PyArrayObject*)PyArray_FromAny(py_ys, 
        PyArray_DescrFromType(NPY_DOUBLE), 1, 1,
        NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED, NULL);
    if (!arr_ys) {
        PyErr_SetString(PyExc_TypeError, "ys must be a 1D contiguous float64 array");
        goto cleanup;
    }

    // 检查数组长度是否一致
    len_x = PyArray_SIZE(arr_xs);
    len_y = PyArray_SIZE(arr_ys);
    if (len_x != len_y) {
        PyErr_SetString(PyExc_ValueError, "xs and ys must have the same length");
        goto cleanup;
    }

    // 检查数组长度是否与 init() 时的 n_srcs 一致，防止 host 侧越界读写
    if (self->obj == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Twinkle object is not initialized");
        goto cleanup;
    }
    if (len_x != (npy_intp)self->obj->n_srcs_all) {
        PyErr_Format(PyExc_ValueError,
            "xs/ys length (%zd) must equal n_srcs in init() (%d)",
            (Py_ssize_t)len_x, self->obj->n_srcs_all);
        goto cleanup;
    }

    try {
        if(RelTol<=0.)      //  如果没设置 RelTol
        {
            // 情况1: ss 是浮点数（标量）
            if (PyFloat_Check(py_ss)) {
                double ss_value = PyFloat_AsDouble(py_ss);
                self->obj->set_params(ss_value, qq, rho, 
                                    (double*)PyArray_DATA(arr_xs), 
                                    (double*)PyArray_DATA(arr_ys));
            }
            // 情况2: ss 是数组
            else {
            arr_ss = (PyArrayObject*)PyArray_FromAny(py_ss, 
                PyArray_DescrFromType(NPY_DOUBLE), 1, 1,
                NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED, NULL);
            
            if (!arr_ss) {
                PyErr_SetString(PyExc_TypeError, 
                    "ss must be float or 1D contiguous float64 array");
                goto cleanup;
            }
            
            // 检查 ss 数组长度
            npy_intp len_s = PyArray_SIZE(arr_ss);
            if (len_s != len_x) {
                PyErr_SetString(PyExc_ValueError, 
                    "ss array length must match xs/ys length");
                goto cleanup;
            }
            
            self->obj->set_params((double*)PyArray_DATA(arr_ss), qq, rho, 
                                 (double*)PyArray_DATA(arr_xs), 
                                 (double*)PyArray_DATA(arr_ys));
            }
        }
        else
        {
            // 情况1: ss 是浮点数（标量）
            if (PyFloat_Check(py_ss)) {
                double ss_value = PyFloat_AsDouble(py_ss);
                self->obj->set_params(ss_value, qq, rho, RelTol,
                                    (double*)PyArray_DATA(arr_xs), 
                                    (double*)PyArray_DATA(arr_ys));
            }
            // 情况2: ss 是数组
            else {
            arr_ss = (PyArrayObject*)PyArray_FromAny(py_ss, 
                PyArray_DescrFromType(NPY_DOUBLE), 1, 1,
                NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED, NULL);
            
            if (!arr_ss) {
                PyErr_SetString(PyExc_TypeError, 
                    "ss must be float or 1D contiguous float64 array");
                goto cleanup;
            }
            
            // 检查 ss 数组长度
            npy_intp len_s = PyArray_SIZE(arr_ss);
            if (len_s != len_x) {
                PyErr_SetString(PyExc_ValueError, 
                    "ss array length must match xs/ys length");
                goto cleanup;
            }
            
            self->obj->set_params((double*)PyArray_DATA(arr_ss), qq, rho, RelTol,
                                 (double*)PyArray_DATA(arr_xs), 
                                 (double*)PyArray_DATA(arr_ys));
            }
        }
        
        // 释放所有数组引用
        Py_XDECREF(arr_xs);
        Py_XDECREF(arr_ys);
        Py_XDECREF(arr_ss);
        Py_RETURN_NONE;
        
    } catch (const std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        // 异常处理也会通过cleanup释放资源
    }

cleanup:
    // 清理所有已分配的数组引用
    Py_XDECREF(arr_xs);
    Py_XDECREF(arr_ys);
    Py_XDECREF(arr_ss);
    return NULL;
}

// return_mag_to()方法包装（需配合NumPy使用）
static PyObject* PyTwinkle_return_mag_to(PyTwinkle *self, PyObject *args) {
    PyArrayObject *py_array;
    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &py_array)) {
        PyErr_SetString(PyExc_TypeError, "Argument must be a numpy array");
        return NULL;
    }

    if (PyArray_NDIM(py_array) != 1 || 
        PyArray_TYPE(py_array) != NPY_DOUBLE) 
    {
        PyErr_SetString(PyExc_ValueError, 
            "Need 1D double array");
        return NULL;
    }

    double* buffer = (double*)PyArray_DATA(py_array);
    const npy_intp size = PyArray_SIZE(py_array);

    if (self->obj == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Twinkle object is not initialized");
        return NULL;
    }
    if (size < (npy_intp)self->obj->n_srcs_all) {
        PyErr_Format(PyExc_ValueError,
            "output array size (%zd) is smaller than n_srcs (%d), would cause out-of-bounds write",
            (Py_ssize_t)size, self->obj->n_srcs_all);
        return NULL;
    }

    try {
        self->obj->return_mag_to(buffer);  // 直接写入NumPy数组
        Py_RETURN_NONE;
    } catch (const std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }
}

// return_Ncross_to()方法包装（需配合NumPy使用）
static PyObject* PyTwinkle_return_Ncross_to(PyTwinkle *self, PyObject *args) {
    PyArrayObject *py_array;
    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &py_array)) {
        PyErr_SetString(PyExc_TypeError, "Argument must be a numpy array");
        return NULL;
    }

    if (PyArray_NDIM(py_array) != 1 || 
        PyArray_TYPE(py_array) != NPY_INT)  // 修改为整型数组
    {
        PyErr_SetString(PyExc_ValueError, 
            "Need 1D integer array");  // 错误信息相应修改
        return NULL;
    }

    int* buffer = (int*)PyArray_DATA(py_array);
    const npy_intp size = PyArray_SIZE(py_array);

    if (self->obj == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Twinkle object is not initialized");
        return NULL;
    }
    if (size < (npy_intp)self->obj->n_srcs_all) {
        PyErr_Format(PyExc_ValueError,
            "output array size (%zd) is smaller than n_srcs (%d), would cause out-of-bounds write",
            (Py_ssize_t)size, self->obj->n_srcs_all);
        return NULL;
    }

    try {
        self->obj->return_Ncross_to(buffer);  // 调用对应的C++函数
        Py_RETURN_NONE;
    } catch (const std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }
}

// return_astrom_to()方法包装（需配合NumPy使用）
static PyObject* PyTwinkle_return_astrom_to(PyTwinkle *self, PyObject *args) {
    PyArrayObject *py_array;
    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &py_array)) {
        PyErr_SetString(PyExc_TypeError, "Argument must be a numpy array");
        return NULL;
    }

    if (PyArray_NDIM(py_array) != 1 || 
        PyArray_TYPE(py_array) != NPY_COMPLEX128)  // 修改为复数数组
    {
        PyErr_SetString(PyExc_ValueError, 
            "Need 1D complex128 array");
        return NULL;
    }

    twinkle::complex_t<double>* buffer = (twinkle::complex_t<double>*)PyArray_DATA(py_array);
    const npy_intp size = PyArray_SIZE(py_array);

    if (self->obj == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Twinkle object is not initialized");
        return NULL;
    }
    if (size < (npy_intp)self->obj->n_srcs_all) {
        PyErr_Format(PyExc_ValueError,
            "output array size (%zd) is smaller than n_srcs (%d), would cause out-of-bounds write",
            (Py_ssize_t)size, self->obj->n_srcs_all);
        return NULL;
    }

    try {
        self->obj->return_astrom_to(buffer);  // 调用对应的C++函数
        Py_RETURN_NONE;
    } catch (const std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }
}

// run()方法包装
static PyObject* PyTwinkle_run(PyTwinkle *self, PyObject *Py_UNUSED(ignored)) {
    try {
        self->obj->run();
        Py_RETURN_NONE;
    } catch (const std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }
}

// run_point_source()方法包装
static PyObject* PyTwinkle_run_pt(PyTwinkle *self, PyObject *Py_UNUSED(ignored)) {
    try {
        self->obj->run_pt();
        Py_RETURN_NONE;
    } catch (const std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }
}

// run_LimbDarkening()方法包装
static PyObject* PyTwinkle_runLD(PyTwinkle *self, PyObject *args) {
    double LD_a;
    
    // 解析参数：期望一个浮点数参数
    if (!PyArg_ParseTuple(args, "d", &LD_a)) {
        return NULL;  // 参数解析失败，Python会自动设置异常
    }
    
    // 可选：验证参数范围
    if (LD_a < 0.0 || LD_a > 1.0) {
        PyErr_SetString(PyExc_ValueError, "LD_a must be between 0.0 and 1.0");
        return NULL;
    }
    
    try {
        self->obj->runLD(LD_a);
        Py_RETURN_NONE;
    } catch (const std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }
}
// static PyObject* PyTwinkle_runLD2(PyTwinkle *self, PyObject *args) {
//     double LD_a;
    
//     // 解析参数：期望一个浮点数参数
//     if (!PyArg_ParseTuple(args, "d", &LD_a)) {
//         return NULL;  // 参数解析失败，Python会自动设置异常
//     }
    
//     // 可选：验证参数范围
//     if (LD_a < 0.0 || LD_a > 1.0) {
//         PyErr_SetString(PyExc_ValueError, "LD_a must be between 0.0 and 1.0");
//         return NULL;
//     }
    
//     try {
//         self->obj->runLD2(LD_a);
//         Py_RETURN_NONE;
//     } catch (const std::exception& e) {
//         PyErr_SetString(PyExc_RuntimeError, e.what());
//         return NULL;
//     }
// }

static PyObject* get_numpy_build_version(PyObject* self, PyObject* args) {
    // 优先使用 NumPy 提供的字符串宏
    #if defined(NPY_VERSION_STRING)
        return PyUnicode_FromString(NPY_VERSION_STRING);
    #else
        // 手动转换版本号
        const int major = NPY_VERSION >> 24;
        const int minor = (NPY_VERSION >> 16) & 0xFF;
        const int patch = NPY_VERSION & 0xFFFF;
        
        char version_str[32];
        snprintf(version_str, sizeof(version_str), "%d.%d.%d", 
                 major, minor, patch);
        
        return PyUnicode_FromString(version_str);
    #endif
}

// dump_margin(source_idx) — copy pool_margin D→H and return boundary data for a single source
static PyObject* PyTwinkle_dump_margin(PyTwinkle *self, PyObject *args) {
    int src_idx = 0;
    if (!PyArg_ParseTuple(args, "i", &src_idx)) {
        PyErr_SetString(PyExc_TypeError, "Usage: dump_margin(source_index)");
        return NULL;
    }
    
    self->obj->dump_margin();
    
    // Find the right source_base for this source index
    int n_src_per_stream = (self->obj->n_srcs_all - 1) / self->obj->vp_sol.size() + 1;
    int stream_idx = src_idx / n_src_per_stream;
    int local_idx  = src_idx % n_src_per_stream;
    
    if (stream_idx >= (int)self->obj->vp_sol.size()) {
        PyErr_SetString(PyExc_IndexError, "source index out of range");
        return NULL;
    }
    
    auto* pool = &self->obj->vp_sol[stream_idx]->pool_margin;
    int n_pts = pool->n_pts;         // n_src * n_point_max
    int n_point_max = pool->stride;  // n_point_max
    int n_src_local = n_pts / n_point_max;
    
    if (local_idx >= n_src_local) {
        PyErr_SetString(PyExc_IndexError, "local index out of range");
        return NULL;
    }
    
    npy_intp npy_npts = (npy_intp)n_point_max;
    
    // Allocate numpy arrays for key fields
    PyObject* arr_Q        = PyArray_SimpleNew(1, &npy_npts, NPY_DOUBLE);
    PyObject* arr_Nphys    = PyArray_SimpleNew(1, &npy_npts, NPY_INT);
    PyObject* arr_skip     = PyArray_SimpleNew(1, &npy_npts, NPY_BOOL);
    PyObject* arr_prev_idx = PyArray_SimpleNew(1, &npy_npts, NPY_INT);
    PyObject* arr_next_idx = PyArray_SimpleNew(1, &npy_npts, NPY_INT);
    PyObject* arr_area_int = PyArray_SimpleNew(1, &npy_npts, NPY_DOUBLE);
    PyObject* arr_err_int  = PyArray_SimpleNew(1, &npy_npts, NPY_DOUBLE);
    PyObject* arr_pos_x    = PyArray_SimpleNew(1, &npy_npts, NPY_DOUBLE);
    PyObject* arr_pos_y    = PyArray_SimpleNew(1, &npy_npts, NPY_DOUBLE);
    
    npy_intp npy_img_dims[2] = {npy_npts, 5};
    PyObject* arr_img_x     = PyArray_SimpleNew(2, npy_img_dims, NPY_DOUBLE);
    PyObject* arr_img_y     = PyArray_SimpleNew(2, npy_img_dims, NPY_DOUBLE);
    PyObject* arr_img_phys  = PyArray_SimpleNew(2, npy_img_dims, NPY_BOOL);
    PyObject* arr_img_parity= PyArray_SimpleNew(2, npy_img_dims, NPY_BOOL);
    PyObject* arr_next_j    = PyArray_SimpleNew(2, npy_img_dims, NPY_INT);
    PyObject* arr_nj_self   = PyArray_SimpleNew(2, npy_img_dims, NPY_INT);
    PyObject* arr_deltaS    = PyArray_SimpleNew(2, npy_img_dims, NPY_DOUBLE);
    PyObject* arr_deltaS_E  = PyArray_SimpleNew(2, npy_img_dims, NPY_DOUBLE);
    PyObject* arr_wedge     = PyArray_SimpleNew(2, npy_img_dims, NPY_DOUBLE);
    
    int base_offset = local_idx * n_point_max;
    
    for (int i = 0; i < n_point_max; i++) {
        auto& pt = pool->dat_h[base_offset + i];
        ((double*)PyArray_DATA(arr_Q))[i]        = (double)pt.Q;
        ((int*)   PyArray_DATA(arr_Nphys))[i]    = pt.Nphys;
        ((bool*)  PyArray_DATA(arr_skip))[i]     = pt.skip;
        ((int*)   PyArray_DATA(arr_prev_idx))[i] = pt.prev_src_idx;
        ((int*)   PyArray_DATA(arr_next_idx))[i] = pt.next_src_idx;
        ((double*)PyArray_DATA(arr_area_int))[i] = (double)pt.area_interval;
        ((double*)PyArray_DATA(arr_err_int))[i]  = (double)pt.error_interval;
        ((double*)PyArray_DATA(arr_pos_x))[i]    = (double)pt.position.re;
        ((double*)PyArray_DATA(arr_pos_y))[i]    = (double)pt.position.im;
        
        for (int j = 0; j < 5; j++) {
            int ij = i * 5 + j;
            ((double*)PyArray_DATA(arr_img_x))[ij]      = (double)pt.images[j].position.re;
            ((double*)PyArray_DATA(arr_img_y))[ij]      = (double)pt.images[j].position.im;
            ((bool*)  PyArray_DATA(arr_img_phys))[ij]   = pt.images[j].physical;
            ((bool*)  PyArray_DATA(arr_img_parity))[ij] = pt.images[j].parity;
            ((int*)   PyArray_DATA(arr_next_j))[ij]     = pt.next_j[j];
            ((int*)   PyArray_DATA(arr_nj_self))[ij]    = pt.next_idx[j];
            ((double*)PyArray_DATA(arr_deltaS))[ij]     = (double)pt.deltaS[j];
            ((double*)PyArray_DATA(arr_deltaS_E))[ij]   = (double)pt.deltaS_Err[j];
            ((double*)PyArray_DATA(arr_wedge))[ij]      = (double)pt.wedge[j];
        }
    }
    
    PyObject* dict = PyDict_New();
    PyDict_SetItemString(dict, "Q",           arr_Q);
    PyDict_SetItemString(dict, "Nphys",       arr_Nphys);
    PyDict_SetItemString(dict, "skip",        arr_skip);
    PyDict_SetItemString(dict, "prev_src_idx",arr_prev_idx);
    PyDict_SetItemString(dict, "next_src_idx",arr_next_idx);
    PyDict_SetItemString(dict, "area_interval",arr_area_int);
    PyDict_SetItemString(dict, "error_interval",arr_err_int);
    PyDict_SetItemString(dict, "pos_x",       arr_pos_x);
    PyDict_SetItemString(dict, "pos_y",       arr_pos_y);
    PyDict_SetItemString(dict, "img_x",       arr_img_x);
    PyDict_SetItemString(dict, "img_y",       arr_img_y);
    PyDict_SetItemString(dict, "img_physical",arr_img_phys);
    PyDict_SetItemString(dict, "img_parity",  arr_img_parity);
    PyDict_SetItemString(dict, "next_j",      arr_next_j);
    PyDict_SetItemString(dict, "next_idx",    arr_nj_self);
    PyDict_SetItemString(dict, "deltaS",      arr_deltaS);
    PyDict_SetItemString(dict, "deltaS_Err",  arr_deltaS_E);
    PyDict_SetItemString(dict, "wedge",       arr_wedge);
    PyDict_SetItemString(dict, "n_point_max", PyLong_FromLong(n_point_max));
    
    Py_DECREF(arr_Q); Py_DECREF(arr_Nphys); Py_DECREF(arr_skip);
    Py_DECREF(arr_prev_idx); Py_DECREF(arr_next_idx);
    Py_DECREF(arr_area_int); Py_DECREF(arr_err_int);
    Py_DECREF(arr_pos_x); Py_DECREF(arr_pos_y);
    Py_DECREF(arr_img_x); Py_DECREF(arr_img_y);
    Py_DECREF(arr_img_phys); Py_DECREF(arr_img_parity);
    Py_DECREF(arr_next_j); Py_DECREF(arr_nj_self);
    Py_DECREF(arr_deltaS); Py_DECREF(arr_deltaS_E); Py_DECREF(arr_wedge);
    
    return dict;
}

// 方法表修正（参数类型匹配）
static PyMethodDef PyTwinkle_methods[] = {
    {"init", (PyCFunction)PyTwinkle_init, METH_VARARGS, "Initialize driver"},
    {"set_params", (PyCFunction)PyTwinkle_set_params, METH_VARARGS, "Set parameters"},
    {"return_mag_to", (PyCFunction)PyTwinkle_return_mag_to, METH_VARARGS, "Write results to numpy array"},
    {"return_Ncross_to", (PyCFunction)PyTwinkle_return_Ncross_to, METH_VARARGS, "Write results to numpy array"},
    {"return_astrom_to", (PyCFunction)PyTwinkle_return_astrom_to, METH_VARARGS, "Write results to numpy array"},
    {"run", (PyCFunction)PyTwinkle_run, METH_NOARGS, "Execute computation"},
    {"run_pt", (PyCFunction)PyTwinkle_run_pt, METH_NOARGS, "point source approximation"},
    {"runLD", (PyCFunction)PyTwinkle_runLD, METH_VARARGS, "Run limb darkening calculation with given linear LD coefficient"},
    {"dump_margin", (PyCFunction)PyTwinkle_dump_margin, METH_VARARGS, "Dump pool_margin D->H for a source index"},
    {"get_numpy_build_version", get_numpy_build_version, METH_NOARGS, "Build-time numpy version"},
    {NULL, NULL, 0, NULL}  // Sentinel
};

// 正确初始化顺序的完整类型定义
static PyTypeObject PyTwinkleType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    /* tp_name */         "twinkle.Twinkle",
    /* tp_basicsize */    sizeof(PyTwinkle),
    /* tp_itemsize */     0,
    /* tp_dealloc */      (destructor)PyTwinkle_dealloc,
    /* tp_vectorcall_offset */ 0,
    /* tp_getattr */      0,
    /* tp_setattr */      0,
    /* tp_as_async */     0,
    /* tp_repr */         0,
    /* tp_as_number */    0,
    /* tp_as_sequence */  0,
    /* tp_as_mapping */   0,
    /* tp_hash */         0,
    /* tp_call */         0,
    /* tp_str */          0,
    /* tp_getattro */     0,
    /* tp_setattro */     0,
    /* tp_as_buffer */    0,
    /* tp_flags */        Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    /* tp_doc */          "Twinkle objects",
    /* tp_traverse */     0,
    /* tp_clear */        0,
    /* tp_richcompare */  0,
    /* tp_weaklistoffset */ 0,
    /* tp_iter */         0,
    /* tp_iternext */     0,
    /* tp_methods */      PyTwinkle_methods,
    /* tp_members */      0,
    /* tp_getset */       0,
    /* tp_base */         0,
    /* tp_dict */         0,
    /* tp_descr_get */    0,
    /* tp_descr_set */    0,
    /* tp_dictoffset */   0,
    /* tp_init */         (initproc)PyTwinkle_init,
    /* tp_alloc */        0,
    /* tp_new */          PyTwinkle_new
};

// 模块定义结构体（新增部分）
static struct PyModuleDef module_def = {
    PyModuleDef_HEAD_INIT,
    .m_name = "twinkle",
    .m_doc = "Twinkle CUDA加速模块",
    .m_size = -1,  // 表示模块状态全局
};

// 模块初始化函数修正版
PyMODINIT_FUNC PyInit_twinkle(void) {
    // 初始化NumPy API（关键！）
    import_array();

    // 创建模块对象
    PyObject* module = PyModule_Create(&module_def);
    if (!module) return NULL;

    // 准备类型对象
    if (PyType_Ready(&PyTwinkleType) < 0) {
        Py_DECREF(module);
        return NULL;
    }

    // 注册类型到模块
    Py_INCREF(&PyTwinkleType);
    if (PyModule_AddObject(module, "Twinkle", (PyObject*)&PyTwinkleType) < 0) {
        Py_DECREF(&PyTwinkleType);
        Py_DECREF(module);
        return NULL;
    }

    return module;
}