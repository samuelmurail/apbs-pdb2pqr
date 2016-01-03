""" Molecule structures for PDB2PQR

    This module contains the base molecular structures for pdb2pqr.

    ----------------------------

    PDB2PQR -- An automated pipeline for the setup, execution, and analysis of Poisson-Boltzmann
    electrostatics calculations

    Copyright (c) 2002-2016, Jens Erik Nielsen, University College Dublin; Nathan A. Baker, Battelle
    Memorial Institute, Developed at the Pacific Northwest National Laboratory, operated by Battelle
    Memorial Institute, Pacific Northwest Division for the U.S. Department Energy; Paul Czodrowski &
    Gerhard Klebe, University of Marburg.

    All rights reserved.

    Redistribution and use in source and binary forms, with or without modification, are permitted
    provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this list of conditions
      and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice, this list of
      conditions and the following disclaimer in the documentation and/or other materials provided
      with the distribution.
    * Neither the names of University College Dublin, Battelle Memorial Institute, Pacific Northwest
      National Laboratory, US Department of Energy, or University of Marburg nor the names of its
      contributors may be used to endorse or promote products derived from this software without
      specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
    IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
    FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
    CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
    OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.

    ----------------------------
"""

__date__ = "2016-01-01"
__author__ = "Todd Dolinsky, Nathan Baker"

from .structures import Residue

class WAT(Residue):
    """ This class gives data about the Water object, and inherits off the base residue class. """
    water_residue_names = ['HOH', 'WAT']
    def __init__(self, atoms, ref):
        sample_atom = atoms[-1]
        self.atoms = []
        self.name = sample_atom.res_name
        self.chain_id = sample_atom.chain_id
        self.res_seq = sample_atom.res_seq
        self.insert_code = sample_atom.insert_code
        self.fixed = 0
        self.ffname = "WAT"
        self.map = {}
        self.reference = ref
        for atom in atoms:
            if atom.name in ref.altnames: # Rename atoms
                atom.name = ref.altnames[atom.name]
            atom_ = Atom(atom, "HETATM", self)
            atomname = atom_.get("name")
            if atomname not in self.map:
                self.add_atom(atom_)
            else: # Don't add duplicate atom with alt_loc field
                oldatom = self.get_atom(atomname)
                oldatom.set("alt_loc", "")
    def create_atom(self, atomname, newcoords):
        """ Create a water atom.  Note the HETATM field.
        Parameters
        * atomname: The name of the atom (string)
        * newcoords:  The new coordinates of the atom (list)"""
        oldatom = self.atoms[0]
        newatom = Atom(oldatom, "HETATM", self)
        newatom.set("x", newcoords[0])
        newatom.set("y", newcoords[1])
        newatom.set("z", newcoords[2])
        newatom.set("name", atomname)
        newatom.set("occupancy", 1.00)
        newatom.set("temperature_factor", 0.00)
        newatom.added = 1
        self.add_atom(newatom)
    def add_atom(self, atom):
        """ Override the existing add_atom - include the link to the reference object """
        self.atoms.append(atom)
        atomname = atom.get("name")
        self.map[atomname] = atom
        try:
            atom.reference = self.reference.map[atomname]
            for bond in atom.reference.bonds:
                if self.has_atom(bond):
                    bondatom = self.map[bond]
                    if bondatom not in atom.bonds:
                        atom.bonds.append(bondatom)
                    if atom not in bondatom.bonds:
                        bondatom.bonds.append(atom)
        except KeyError:
            atom.reference = None

class LIG(Residue):
    """ This class gives data about the generic ligand object, and inherits off the base residue
    class. """
    def __init__(self, atoms, ref):
        sample_atom = atoms[-1]
        self.atoms = []
        self.name = sample_atom.res_name
        self.chain_id = sample_atom.chain_id
        self.res_seq = sample_atom.res_seq
        self.insert_code = sample_atom.insert_code
        self.fixed = 0
        self.ffname = "WAT"
        self.map = {}
        self.reference = ref
        self.is_n_terminus = False
        self.is_c_terminus = False
        for atom in atoms:
            if atom.name in ref.altnames: # Rename atoms
                atom.name = ref.altnames[atom.name]
            atom_ = Atom(atom, "HETATM", self)
            atomname = atom_.get("name")
            if atomname not in self.map:
                self.add_atom(atom_)
            else: # Don't add duplicate atom with alt_loc field
                oldatom = self.get_atom(atomname)
                oldatom.set("alt_loc", "")
    def create_atom(self, atomname, newcoords):
        """ Create a water atom.  Note the HETATM field.
        Parameters
        * atomname: The name of the atom (string)
        * newcoords:  The new coordinates of the atom (list) """
        oldatom = self.atoms[0]
        newatom = Atom(oldatom, "HETATM", self)
        newatom.set("x", newcoords[0])
        newatom.set("y", newcoords[1])
        newatom.set("z", newcoords[2])
        newatom.set("name", atomname)
        newatom.set("occupancy", 1.00)
        newatom.set("temperature_factor", 0.00)
        newatom.added = 1
        self.add_atom(newatom)
    def add_atom(self, atom):
        """ Override the existing add_atom - include the link to the reference object """
        self.atoms.append(atom)
        atomname = atom.get("name")
        self.map[atomname] = atom
        try:
            atom.reference = self.reference.map[atomname]
            for bond in atom.reference.bonds:
                if self.has_atom(bond):
                    bondatom = self.map[bond]
                    if bondatom not in atom.bonds:
                        atom.bonds.append(bondatom)
                    if atom not in bondatom.bonds:
                        bondatom.bonds.append(atom)
        except KeyError:
            atom.reference = None
