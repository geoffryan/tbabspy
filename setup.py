import setuptools
import subprocess
from pathlib import Path
from setuptools import setup, Extension
# from numpy.distutils.core import setup as npsetup
# from numpy.distutils.core import Extension as npExtension
# from numpy.distutils.
from setuptools.command.build_ext import build_ext
from distutils.errors import DistutilsSetupError
from distutils import log as distutils_logger

import numpy as np

version = {}
with open("tbabspy/version.py", "r") as f:
    exec(f.read(), version)

with open("README.md", "r") as f:
    long_description = f.read()

coreInc = [np.get_include()]

coreSrc = ["tbabspy/tbabscoremodule.c", "tbabspy/util.c",
           # "tbabspy/phfit2Src/phfit2.f",
           "tbabspy/tbabsSrc/tbvabs.c"]
coreDep = ["tbabspy/util.h",
           "tbabspy/tbabsSrc/tbvabs.h",
           "tbabspy/tbabsSrc/fephoto.h",
           "tbabspy/tbabsSrc/nephoto.h",
           "tbabspy/tbabsSrc/ophoto.h"]

phfit2Src = "tbabspy/phfit2Src/phfit2.f"
phfit2Lib = "tbabspy/libphfit2.so"

coreLib = ['phfit2']



# phfit2Src = ["tbabspy/phfit2Src/phfit2.f"]

# phfit2Module = Extension('tbabspy.phfit2', sources=phfit2Src,
                         # extra_compile_args=['-Wall'])

coreModule = Extension('tbabspy.tbabscore', sources=coreSrc,
                        depends=coreDep, libraries=coreLib,
                        library_dirs=['tbabspy'],
                        include_dirs=coreInc,
                        extra_compile_args=['-Wall'])

# grabbed from lucianopaz's answer on
# https://stackoverflow.com/questions/41169711/
#  python-setuptools-distutils-custom-build-for-the-extra-package-with-makefile
class custom_build_ext(build_ext, object):

    my_ext = coreModule.name

    def build_extension(self, ext):

        # If you're not building tbabscore, do the usual thing
        if ext.name != self.my_ext:
            super(custom_build_ext, self).build_extension(ext)

        else:
            """
            sources = ext.sources
            if sources is None or not isinstance(sources, (list, tuple)):
                raise DistutilsSetupError(
                        "in 'ext_modules' option (extension '%s'), "
                        "'sources' must be present and must be "
                        "a list of source filenames" % ext.name)
            sources = list(sources)
            """
            print("OH HOH HOH")

            distutils_logger.info("Compiling " + phfit2Lib)

            subprocess.call(["gfortran", "-O3", "-shared", "-fpic", "-Wall",
                             "-o", phfit2Lib,
                             "-c", phfit2Src])
            
            distutils_logger.info("Compiled " + phfit2Lib)

            super(custom_build_ext, self).build_extension(ext)



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
    ext_modules=[coreModule],
    cmdclass={'build_ext': custom_build_ext},
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
