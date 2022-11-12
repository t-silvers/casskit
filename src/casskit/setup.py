from setuptools import setup, find_namespace_packages

setup(
    name='io_tcga',
    packages=find_namespace_packages(include=['io.*']),
    zip_safe=False,
)