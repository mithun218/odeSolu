# coding: utf-8
import os
import fnmatch
import sysconfig

from setuptools import setup, find_packages
from setuptools.command.build_py import build_py as _build_py

from Cython.Build import cythonize


EXCLUDE_FILES = [
    'odeSolu/seriesSolu.py'
]


# noinspection PyPep8Naming
class build_py(_build_py):

    def find_package_modules(self, package, package_dir):
        ext_suffix = sysconfig.get_config_var('EXT_SUFFIX')
        modules = super().find_package_modules(package, package_dir)
        filtered_modules = []
        for (pkg, mod, filepath) in modules:
            if os.path.exists(filepath.replace('.pyx', ext_suffix)) and os.path.exists(filepath.replace('.py', ext_suffix)):
                continue
            filtered_modules.append((pkg, mod, filepath, ))
        return filtered_modules


def get_ext_paths(root_dir, exclude_files):
    """get filepaths for compilation"""
    paths = []

    for root, dirs, files in os.walk(root_dir):
        for filename in files:
            if os.path.splitext(filename)[1] != '.pyx':
                continue

            file_path = os.path.join(root, filename)
            if file_path in exclude_files:
                continue

            paths.append(file_path)
    return paths


setup(
    name='odeSolu',
    version='1.0.0',
    description='A Python package to get power series solution of a nonlinear ODE ',
    
    packages=find_packages(),
    ext_modules=cythonize(
        get_ext_paths('odeSolu', EXCLUDE_FILES),
        compiler_directives={'language_level': 3}
    ),

    author='Mithun Bairagi',
    author_email='bairagirasulpur@gmail.com',

    license='MIT',

    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        ],
    python_requires='>=3.8',
    install_requires=[
        'numpy>=1.14.5',
        'sympy>=1.0.0'
    ],
    cmdclass={
        'build_py': build_py
    }
)
