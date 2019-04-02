import os

from setuptools import find_packages, setup

here = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="connectf-backend",
    description='connectf backend',
    version='1.0.0',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Zachary Juang',
    classifiers=[
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7'
    ],
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'pandas',
        'numpy',
        'scipy',
        'django',
        'matplotlib',
        'pyparsing',
        'seaborn',
        'lxml',
        'networkx',
        'django-cors-headers',
        'XlsxWriter',
        'mysqlclient',
        'statsmodels',
        'scikit-learn',
        'gunicorn',
        'pyyaml'
    ],
    python_requires='>=3.5'
)
