from setuptools import setup, Extension
import numpy as np

version = {}
with open("tbabspy/version.py", "r") as f:
    exec(f.read(), version)

with open("README.md", "r") as f:
    long_description = f.read()

inc = [np.get_include()]
libs = []
libdirs = []

coreSrc = ["tbabspy/tbabscoremodule.c", "tbabspy/util.c",
           "tbabspy/phfit2Src/phfit2.f",
           "tbabspy/tbabsSrc/tbvabs.c"]
coreDep = ["tbabspy/util.h",
           "tbabspy/tbabsSrc/tbvabs.h",
           "tbabspy/tbabsSrc/fephoto.h",
           "tbabspy/tbabsSrc/nephoto.h",
           "tbabspy/tbabsSrc/ophoto.h"]

tbabscoreModule = Extension('tbabspy.tbabscore', sources=coreSrc,
                            depends=coreDep,
                            extra_compile_args=['-Wall'])

setup(
    name='tbabspy',
    version=version['__version__'],
    author="Geoffrey Ryan",
    author_email="gsryan@umd.edu",
    description='tbabs X-ray absorption models',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/geoffryan/tbabspy',
    packages=['tbabspy'],
    ext_modules=[tbabscoreModule],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: C",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy"],
    install_requires=['numpy>=1.10'],
    extras_require={
        'docs': ['numpydoc']
        },
    project_urls={
        "Source Code": "https://github.com/geoffryan/tbabspy",
        "tbabs XSPEC manual": ("https://heasarc.gsfc.nasa.gov/xanadu/"
                               +"xspec/manual/node265.html"),
        "tbabs beta development": ("https://pulsar.sternwarte.uni-erlangen.de"
                                   +"/wilms/research/tbabs/"),
        "phfit2 Code": "http://www.pa.uky.edu/~verner/photo.html"}
    )
