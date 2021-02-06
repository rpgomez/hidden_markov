from setuptools import setup, find_packages

setup(name="hmm",version="0.2.1", packages=find_packages(),
            py_modules = ['hmm'],
      install_requires = ['numpy','scipy','tqdm','numba'])
