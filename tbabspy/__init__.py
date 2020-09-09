#!/usr/bin/env python3
"""
=======
tbabspy
=======

X-ray absorption in Python
"""
from .version import __version__
from . import tbabscore
from .tbabs import tbabs, dgami

__all__ = ['__version__',
           'tbabscore',
           'tbabs', 'dgami']

