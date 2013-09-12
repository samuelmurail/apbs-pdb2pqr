###############################################################################
#
# build_utils.py
#
# Builders and utilities to support scons building.
#
###############################################################################

from SCons.Script import *
import SCons.Errors
import re

def CompilePythonAction(targetfile, sourcefile):
    """Compile python into byte code.
    """
    try:
        import py_compile
    except ImportError:
        raise SCons.Errors.InternalError, "Couldn't import py_compile module"
    
    try:
        py_compile.compile(sourcefile, targetfile, doraise=True)
    except py_compile.PyCompileError:
        raise SCons.Errors.InternalError, "Couldn't compile {0}".format(sourcefile)
    
    

def CompilePythonActionStringFunc(targetfile, sourcefile):
    return 'Compiling python to bytecode ("%s", "%s")' % (targetfile, sourcefile)

CompilePython = SCons.Action.ActionFactory( CompilePythonAction, CompilePythonActionStringFunc )




