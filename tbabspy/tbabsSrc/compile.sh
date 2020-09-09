#!/bin/sh

echo "initpackage absmodel lmodel_tbnew.dat . \nquit\ny" | xspec
echo "lmod absmodel .\nquit\ny" | xspec
rm -f *~ *.o *FunctionMap.* lpack_* Makefile

