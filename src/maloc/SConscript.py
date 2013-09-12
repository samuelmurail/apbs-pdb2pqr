import fnmatch
import os
Import('env')

CPPPATH_ext = []

for root, dirnames, filenames in os.walk('src'):
  for dirname in fnmatch.filter(dirnames, 'maloc'):
      env.Append(CPPPATH=[os.path.abspath(root)])

sourceObjects=[]
for root, dirnames, filenames in os.walk('src'):
  for filename in fnmatch.filter(filenames, '*.c'):
      o = env.Object(os.path.join(root, filename), CPPPATH=[root, 'src']+env['CPPPATH'])
      sourceObjects.append(o)


staticLib = env.StaticLibrary(target='../../lib/maloc', source=sourceObjects)

Return('staticLib')
