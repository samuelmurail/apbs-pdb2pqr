from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import os


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


class build_py(build_py):
    def byte_compile(self, files):
        super(build_py, self).byte_compile(files)
        for file in files:
            if file.endswith('.py'):
                tmp = os.path.split(file)
                if tmp[0] == "build/lib/pdb2pqr" and not tmp[1] in ['__init__.py', 'main.py']:
                    print(file + "<--- DELETED")
                    os.unlink(file)
                if tmp[0] == "build/lib/pdb2pqr/pdb2pka":
                    print(file + "<--- DELETED")
                    os.unlink(file)


setup(name='pdb2pqr',
      version='0.0.1',
      url='http://www.poissonboltzmann.org/',
      description="PDB2PQR: an automated pipeline for the setup of Poisson-Boltzmann electrostatics calculations",
      long_description=read("pdb2pqr/README.md"),
      license="BSD",
      packages=['pdb2pqr',
                'pdb2pqr.src',
                'pdb2pqr.pdb2pka',
                'pdb2pqr.propka30',
                'pdb2pqr.propka30.Source',
                'pdb2pqr.extensions'
                ],
      package_data={'pdb2pqr': ['dat/*', 'NEWS', 'README.md', 'COPYING', 'AUTHORS', 'propka30/Source/*.dat', 'propka30/Source/ions.list']},
      cmdclass=dict(build_py=build_py),
      entry_points={
          'console_scripts': [
              'pdb2pqr_cli = pdb2pqr.main:main'
          ]
      }
      )
