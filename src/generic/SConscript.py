import fnmatch
import os
Import('env')

sourceObjects=[]
for root, dirnames, filenames in os.walk('.'):
  for filename in fnmatch.filter(filenames, '*.c'):
      o = env.Object(os.path.join(root, filename))
      sourceObjects.append(o)


staticLib = env.StaticLibrary(target='build/apbs', source=sourceObjects)

Return('staticLib')
