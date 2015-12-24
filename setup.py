from setuptools import setup, find_packages

from version import get_git_version


install_requires = ['numpy', 'matplotlib']

setup(
    name='ngs-tk',
    version=get_git_version(),
    url=None,
    author='Julian de Ruiter',
    author_email='julianderuiter@gmail.com',
    description='Toolkit for various NGS analyses.',
    license='BSD',
    packages=find_packages(),
    include_package_data=True,
    zip_safe=True,
    classifiers=[],
    install_requires=install_requires
)
