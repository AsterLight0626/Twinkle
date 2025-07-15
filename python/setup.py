# python setup.py build_ext --inplace

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import os
import sys
import numpy as np
from sysconfig import get_paths as sysconfig_get_paths

python_include = sysconfig_get_paths()['include']
CUDA_HOME = os.environ.get("CUDA_HOME", "/usr/local/cuda")

class CUDAExtension(Extension):
    def __init__(self, name, sources, **kwargs):
        kwargs.setdefault('extra_compile_args', [])
        kwargs.setdefault('include_dirs', [])
        super().__init__(name, sources, **kwargs)

class CustomBuildExt(build_ext):
    def build_extension(self, ext):
        if isinstance(ext, CUDAExtension):
            objects = self._compile_cuda(ext)
            
            # 获取输出路径并创建目录
            ext_path = self.get_ext_fullpath(ext.name)
            ext_dir = os.path.dirname(ext_path)
            os.makedirs(ext_dir, exist_ok=True)

            # 使用NVCC进行最终链接（关键修改）
            link_cmd = [
                '/usr/local/cuda-12.3/bin/nvcc',
                '--shared',
                '-arch=sm_70',
                '-Xcompiler', '-fPIC',       # 显式传递PIC参数
                '-Xlinker', '-rpath=/usr/local/cuda-12.3/lib64',
                '-o', ext_path,
                *objects,
                '-L/usr/local/cuda-12.3/lib64',
                '-lcudart',
                '-lcudadevrt',
                '-Xcompiler', '-lstdc++'     # 传递C++标准库
            ]
            
            print("Final link command:", ' '.join(link_cmd))
            self.spawn(link_cmd)
        else:
            super().build_extension(ext)

    def _compile_cuda(self, ext):
        nvcc = os.path.join(CUDA_HOME, "bin", "nvcc")
        if not os.path.exists(nvcc):
            raise RuntimeError(f"NVCC not found at {nvcc}")

        nvcc_flags = [
            '-arch=sm_70',
            '-std=c++17',
            '-O3',
            '-x', 'cu',
            '-rdc=true',
            '--expt-relaxed-constexpr',
            '-Xcompiler', '-fPIC',          # 主机代码PIC
            '-Xcompiler', '-D_GLIBCXX_USE_CXX11_ABI=0',
            '-Xcompiler', '-DNPY_NO_DEPRECATED_API=NPY_2_0_API_VERSION',  # 新增这行
            '-Xcompiler', '-Wno-deprecated-declarations'  # 抑制警告
        ]

        include_dirs = [
            os.path.abspath('../src'),
            np.get_include(),
            os.path.join(CUDA_HOME, 'include'),
            python_include
        ]

        objects = []
        for source in ext.sources:
            obj_path = os.path.join(
                self.build_temp,
                os.path.splitext(os.path.basename(source))[0] + ".o"
            )
            os.makedirs(os.path.dirname(obj_path), exist_ok=True)
            
            compile_cmd = [
                nvcc,
                '-c', source,
                '-o', obj_path,
                *[f'-I{d}' for d in include_dirs],
                *nvcc_flags
            ]
            self.spawn(compile_cmd)
            objects.append(obj_path)
        
        return objects

src_files = []
for root, dirs, files in os.walk(os.path.join('..', 'src')):
    for file in files:
        if file.endswith(('.cu', '.cpp')) and not file.endswith("main.cpp"):
            src_files.append(os.path.join(root, file))

cuda_extension = CUDAExtension(
    name='twinkle',
    sources=['wrapper.cpp'] + src_files,
    include_dirs=[
        '../src',
        np.get_include(),
        os.path.join(CUDA_HOME, 'include'),
        python_include
    ],
    library_dirs=['/usr/local/cuda-12.3/lib64'],
    language='c++'
)
# print("NumPy include path:", np.get_include())

setup(
    name='twinkle',
    version='1.0',
    install_requires=['numpy'],
    ext_modules=[cuda_extension],
    cmdclass={'build_ext': CustomBuildExt}
)