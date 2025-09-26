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

    // 定义参数选项
    static char* kwlist1[] = {(char*)"n_srcs", (char*)"device_num", (char*)"n_stream", (char*)"RelTol", NULL};
    static char* kwlist2[] = {(char*)"n_srcs", (char*)"device_num", (char*)"n_stream", NULL};
    
    // 解析参数字典
    if (PyArg_ParseTupleAndKeywords(args, kwargs, "iiid", kwlist1, &n_srcs, &device_num, &n_stream, &RelTol))
    {    }
    else
    {
        PyErr_Clear();
        if (PyArg_ParseTupleAndKeywords(args, kwargs, "iii", kwlist2, &n_srcs, &device_num, &n_stream))
        {  
            printf("In new version, RelTol is set in init(), default 1e-4. The RelTol in set_params() won't be used\n");
        }
        else
        {
            return -1;
        }
    }


    try {
        self->obj = new twinkle::driver_t();
        self->obj->init(n_srcs, device_num, n_stream, RelTol);
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

// // set_params()方法包装
// static PyObject* PyTwinkle_set_params(PyTwinkle *self, PyObject *args) {
//     double ss, qq, rho, RELTOL;
//     PyObject *py_xs, *py_ys;

//     // 解析参数：ss, qq, rho, RELTOL, xs(list), ys(list)
//     if (!PyArg_ParseTuple(args, "ddddO!O!", 
//                          &ss, &qq, &rho, &RELTOL,
//                          &PyList_Type, &py_xs,
//                          &PyList_Type, &py_ys)) 
//     {
//         PyErr_SetString(PyExc_TypeError, 
//             "Arguments: (ss, qq, rho, RELTOL: float, xs: list, ys: list)");
//         return NULL;
//     }

//     // 转换Python列表到C数组
//     const Py_ssize_t len_x = PyList_Size(py_xs);
//     const Py_ssize_t len_y = PyList_Size(py_ys);
//     if (len_x != len_y) {
//         PyErr_SetString(PyExc_ValueError, "xs and ys must have same length");
//         return NULL;
//     }

//     std::vector<double> xs(len_x), ys(len_y);
//     for (Py_ssize_t i = 0; i < len_x; ++i) {
//         xs[i] = PyFloat_AsDouble(PyList_GetItem(py_xs, i));
//         ys[i] = PyFloat_AsDouble(PyList_GetItem(py_ys, i));
//     }

//     try {
//         self->obj->set_params(ss, qq, rho, RELTOL, xs.data(), ys.data());
//         Py_RETURN_NONE;
//     } catch (const std::exception& e) {
//         PyErr_SetString(PyExc_RuntimeError, e.what());
//         return NULL;
//     }
// }

// static PyObject* PyTwinkle_set_params(PyTwinkle *self, PyObject *args) {
//     double ss, qq, rho, RELTOL;
//     PyObject *py_xs, *py_ys;

//     // 解析参数：ss, qq, rho, RELTOL, xs(array), ys(array)
//     if (!PyArg_ParseTuple(args, "ddddOO", 
//                          &ss, &qq, &rho, &RELTOL,
//                          &py_xs, &py_ys)) {
//         PyErr_SetString(PyExc_TypeError, 
//             "Arguments: (ss, qq, rho, RELTOL: float, xs: np.ndarray, ys: np.ndarray)");
//         return NULL;
//     }

//     PyArrayObject *arr_xs = NULL;
//     PyArrayObject *arr_ys = NULL;

//     // 转换并检查 xs 数组：1维，双精度，C连续
//     arr_xs = (PyArrayObject*) PyArray_FromAny(py_xs, 
//         PyArray_DescrFromType(NPY_DOUBLE), 1, 1,
//         NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED, NULL);
//     if (!arr_xs) {
//         PyErr_SetString(PyExc_TypeError, "xs must be a 1D contiguous float64 array");
//         Py_XDECREF(arr_xs);
//         Py_XDECREF(arr_ys);
//         return NULL;
//     }

//     // 转换并检查 ys 数组
//     arr_ys = (PyArrayObject*) PyArray_FromAny(py_ys, 
//         PyArray_DescrFromType(NPY_DOUBLE), 1, 1,
//         NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED, NULL);
//     if (!arr_ys) {
//         PyErr_SetString(PyExc_TypeError, "ys must be a 1D contiguous float64 array");
//         Py_XDECREF(arr_xs);
//         Py_XDECREF(arr_ys);
//         return NULL;
//     }

//     // 检查数组长度是否一致
//     npy_intp len_x = PyArray_SIZE(arr_xs);
//     npy_intp len_y = PyArray_SIZE(arr_ys);
//     if (len_x != len_y) {
//         PyErr_SetString(PyExc_ValueError, "xs and ys must have the same length");
//         Py_XDECREF(arr_xs);
//         Py_XDECREF(arr_ys);
//         return NULL;
//     }

//     // 获取数据指针
//     double *xs_data = (double*) PyArray_DATA(arr_xs);
//     double *ys_data = (double*) PyArray_DATA(arr_ys);

//     try {
//         // 调用 C++ 函数，传递原始指针
//         self->obj->set_params(ss, qq, rho, RELTOL, xs_data, ys_data);
//         Py_DECREF(arr_xs);
//         Py_DECREF(arr_ys);
//         Py_RETURN_NONE;
//     } catch (const std::exception& e) {
//         PyErr_SetString(PyExc_RuntimeError, e.what());
//         Py_XDECREF(arr_xs);
//         Py_XDECREF(arr_ys);
//         return NULL;
//     }

// }

// static PyObject* PyTwinkle_set_params(PyTwinkle *self, PyObject *args) {
//     PyObject *py_ss;  // 通用对象，可以是浮点数或数组
//     double qq, rho, RelTol;
//     PyObject *py_xs, *py_ys;

//     // 解析参数：ss(标量或数组), qq, rho, RELTOL, xs(array), ys(array)
//     if (PyArg_ParseTuple(args, "OdddOO", &py_ss, &qq, &rho, &RelTol, &py_xs, &py_ys)) {
//         // 第一种情况：提供了RelTol
//         // printf("In new version, RelTol is set in Twinkle.init(), default 1e-4. The RelTol here won't be used\n");
//     }
//     else {
//         PyErr_Clear();
//         // 尝试第二种参数组合（不包含RelTol）
//         if (!PyArg_ParseTuple(args, "OddOO", &py_ss, &qq, &rho, &py_xs, &py_ys)) {
//             PyErr_SetString(PyExc_TypeError, 
//                 "Arguments: (ss: float OR np.ndarray, qq, rho: float, "
//                 "xs: np.ndarray, ys: np.ndarray) OR "
//                 "(ss: float OR np.ndarray, qq, rho: float, RelTol: float, "
//                 "xs: np.ndarray, ys: np.ndarray)");
//             return NULL;
//         }
//         // 第二种情况：不提供RelTol
//     }

//     // 转换并检查 xs 数组
//     PyArrayObject *arr_xs = (PyArrayObject*)PyArray_FromAny(py_xs, 
//         PyArray_DescrFromType(NPY_DOUBLE), 1, 1,
//         NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED, NULL);
//     if (!arr_xs) {
//         PyErr_SetString(PyExc_TypeError, "xs must be a 1D contiguous float64 array");
//         return NULL;
//     }

//     // 转换并检查 ys 数组
//     PyArrayObject *arr_ys = (PyArrayObject*)PyArray_FromAny(py_ys, 
//         PyArray_DescrFromType(NPY_DOUBLE), 1, 1,
//         NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED, NULL);
//     if (!arr_ys) {
//         PyErr_SetString(PyExc_TypeError, "ys must be a 1D contiguous float64 array");
//         return NULL;
//     }

//     // 检查数组长度是否一致
//     npy_intp len_x = PyArray_SIZE(arr_xs);
//     npy_intp len_y = PyArray_SIZE(arr_ys);
//     if (len_x != len_y) {
//         PyErr_SetString(PyExc_ValueError, "xs and ys must have the same length");
//         return NULL;
//     }

//     // 获取数据指针
//     double *xs_data = (double*)PyArray_DATA(arr_xs);
//     double *ys_data = (double*)PyArray_DATA(arr_ys);
    
//     try {
//         // 情况1: ss 是浮点数（标量）
//         if (PyFloat_Check(py_ss)) {
//             double ss_value = PyFloat_AsDouble(py_ss);
//             self->obj->set_params(ss_value, qq, rho, xs_data, ys_data);
//         }
//         // 情况2: ss 是数组
//         else {
//             PyArrayObject *arr_ss = (PyArrayObject*)PyArray_FromAny(py_ss, 
//                 PyArray_DescrFromType(NPY_DOUBLE), 1, 1,
//                 NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED, NULL);
            
//             if (!arr_ss) {
//                 PyErr_SetString(PyExc_TypeError, 
//                     "ss must be float or 1D contiguous float64 array");
//                 return NULL;
//             }
            
//             // 检查 ss 数组长度
//             npy_intp len_s = PyArray_SIZE(arr_ss);
//             if (len_s != len_x) {
//                 PyErr_SetString(PyExc_ValueError, 
//                     "ss array length must match xs/ys length");
//                 return NULL;
//             }
            
//             double *ss_data = (double*)PyArray_DATA(arr_ss);
//             self->obj->set_params(ss_data, qq, rho, xs_data, ys_data);
//         }
        
//         Py_RETURN_NONE;
        
//     } catch (const std::exception& e) {
//         PyErr_SetString(PyExc_RuntimeError, e.what());
//         return NULL;
//     }
// }

static PyObject* PyTwinkle_set_params(PyTwinkle *self, PyObject *args) {
    PyObject *py_ss = nullptr;
    double qq = 0.0, rho = 0.0, RelTol = 0.0;
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

    try {
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

    try {
        self->obj->return_mag_to(buffer);  // 直接写入NumPy数组
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

// 方法表修正（参数类型匹配）
static PyMethodDef PyTwinkle_methods[] = {
    {"init", (PyCFunction)PyTwinkle_init, METH_VARARGS, "Initialize driver"},
    {"set_params", (PyCFunction)PyTwinkle_set_params, METH_VARARGS, "Set parameters"},
    {"return_mag_to", (PyCFunction)PyTwinkle_return_mag_to, METH_VARARGS, "Write results to numpy array"},
    {"run", (PyCFunction)PyTwinkle_run, METH_NOARGS, "Execute computation"},
    {"runLD", (PyCFunction)PyTwinkle_runLD, METH_VARARGS, "Run limb darkening calculation with given linear LD coefficient"},
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