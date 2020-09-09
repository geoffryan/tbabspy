#!/usr/bin/env python3
"""
=======
tbabspy
=======

X-ray absorption in Python
"""
from .version import __version__
from . import tbabscore
from . import util
from .tbabs import tbabs, ztbabs, tbfeo, tbgas

__all__ = ['__version__',
           'tbabscore', 'util',
           'tbabs', 'ztbabs', 'tbfeo', 'tbgas',
           'dgami', 'phfit2', 'vernabs']

