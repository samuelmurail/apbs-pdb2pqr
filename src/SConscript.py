import fnmatch
import os
Import('env')

libs = [SConscript('maloc/SConscript.py')]

def build_stuff(dirname):
    sourceObjects=[]
    for root, dirnames, filenames in os.walk(dirname):
      for filename in fnmatch.filter(filenames, '*.c'):
         o = env.Object(os.path.join(root, filename))
         sourceObjects.append(o)


    return env.StaticLibrary(target='../lib/apbs_'+dirname, source=sourceObjects)

libs.append(build_stuff('generic'))
libs.append(build_stuff('mg'))
libs.append(build_stuff('pmgc'))

#Cus F you circular dependencies.
apbs = env.Program(target='../bin/apbs', source=['main.c', 'routines.c'], LIBS=(libs*2)+env['LIBS'])

Default(apbs)
