import distutils.sysconfig

import os
from defaults import *
import atexit

config_file = 'build_config.py'

gcv = distutils.sysconfig.get_config_var

vars = Variables(['.variables.cache', config_file], ARGUMENTS)

vars.Add(PathVariable('PREFIX',
                      'Install directory',
                      defaultPrefix,
                      PathVariable.PathAccept))
					  
vars.Add(BoolVariable('REBUILD_SWIG', 
					  'Set to True to rebuild the swig bindings. Requires swig on the the user path.',
					  False))
if os.name == 'nt':
    tool_chain = ['mingw']
else:
    tool_chain = ['default']
    
env = Environment(variables=vars,
                  tools=tool_chain + ['swig', 'textfile'], 
                  SWIGFLAGS=['-python', '-c++'], 
                  SHLIBPREFIX="", 
                  SHLIBSUFFIX=gcv('SO'),
                  LDMODULESUFFIX=gcv('SO'),
                  LIBS=[])


#python_lib = 'python' + gcv('VERSION')
#env.Append(LIBS=[python_lib])
#env.Append(ENV={'PATH' : os.environ['PATH']})
env.Append(CPPPATH=[os.path.abspath('src')])
env.Append(CFLAGS=['-std=c99'])

#if os.name == 'nt':
#    python_root = sys.prefix    
#    python_include = os.path.join(python_root, 'include')
#    python_libs = os.path.join(python_root, 'libs')    
#    env.Append(LIBPATH=[python_libs])
#else:
#    env.Append(LIBPATH=[gcv('LIBDIR')])
    
Export('env')

prefix = env['PREFIX']
prefix = prefix.replace('\\', '/')
if prefix[-1] != '/':
    prefix+='/'
    
env['PREFIX'] = prefix

Help(vars.GenerateHelpText(env))

vars.Save('.variables.cache', env)

#Not the optimal way to do this...
#Should figure out how to do it with a delete command
Clean('bin/apbs', '.variables.cache')

if not GetOption("clean"):
    SConscript('SConscript-config.py')

env.Clean("distclean",
          [
           ".sconsign.dblite",
           ".sconf_temp",
           "config.log",
          ])

SConscript('src/SConscript.py')

SConscript('SConscript-error.py')


def print_default_message(target_list):
    target_list = map(str, target_list)
    if any('test' in x for x in target_list):
        return
    if GetOption("clean"):
        return

    print "TODO: Add postbuild help"
    
    
atexit.register(print_default_message, BUILD_TARGETS)
