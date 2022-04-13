#!/usr/bin/python

# Algorithm for Tracking Convective Systems (ATRACKCS)

from setuptools import find_packages, setup
#from distutils.core import setup

setup(name='atrackcs',
        version='1.0',
        description='ATRACKCS is a Python package for the detection and tracking of convective systems using image-processing techniques.',
        author='Álvaro Ramírez Cardona',
        author_email='alramirezca@unal.edu.co',
        url='https://github.com/alramirezca/atrackcs',
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Natural Language :: English',
            'Operating System :: POSIX :: Linux',
            'Operating System :: MacOS',
            'Operating System :: Microsoft :: Windows',
            'Programming Language :: Python :: 3',
            'Topic :: Scientific/Engineering :: Atmospheric Science'
            ],
        install_requires=[
            "bokeh==2.3.2",
            "bottleneck==1.3.2",
            "cartopy==0.19.0.post1",
            "certifi",
            "cftime==1.4.1",
            "click==7.1.2",
            "cloudpickle==1.6.0",
            "cytoolz==0.11.0", 
	    "dask==2021.4.1", 
	    "fsspec==2021.4.0", 
	    "gdal==3.1.4", 
	    "jinja2==2.11.3", 
	    "locket==0.2.0",
	    "markupsafe==1.1.1", 
	    "matplotlib==3.3.4",
	    "mkl-service==2.3.0", 
	    "netcdf4==1.5.6",
	    "numpy==1.20.1",
            "packaging==20.9",
	    "pillow==8.2.0",
	    "pip==21.0.1",
	    "python-dateutil==2.8.1",
	    "sip==4.19.13",
     	    "six==1.15.0",
            "sortedcontainers==2.3.0",
            "typing-extensions==3.7.4.3",
	    "imageio", 
	    "rasterio",
	    "geopandas",
	    "rtree==0.9.3",
  	    "rioxarray",
	    "maup",
	    "folium", 
            "geopy"
        ],
        packages=find_packages(include=['atrackcs', 'atrackcs.*']),
        license='GPL-3.0-or-later'
        )