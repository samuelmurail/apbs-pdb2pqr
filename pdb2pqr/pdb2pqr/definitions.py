""" This file contains classes associated with Amino Acid and Rotamer definitions as used by
PDB2PQR.

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

* Redistributions of source code must retain the above copyright notice, this list of conditions and
  the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, this list of conditions
  and the following disclaimer in the documentation and/or other materials provided with the
  distribution.
* Neither the names of University College Dublin, Battelle Memorial Institute, Pacific Northwest
  National Laboratory, US Department of Energy, or University of Marburg nor the names of its
  contributors may be used to endorse or promote products derived from this software without
  specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

----------------------------
"""

__date__ = "2016-01-02"
__author__ = "Jens Erik Nielsen, Todd Dolinsky, Yong Huang"

AAPATH = "dat/AA.xml"
NAPATH = "dat/NA.xml"
PATCHPATH = "dat/PATCHES.xml"

import copy
import re
from xml import sax
from .utilities import get_data_file
from .structures import Residue, Atom, BACKBONE

from .errors import PDBInternalError

class DefinitionHandler(sax.ContentHandler):
    """ Handler for definitions """
    def __init__(self):
        super().__init__()
        self.curelement = ""
        self.curatom = None
        self.curholder = None
        self.curobj = None
        self.map = {}
        self.patches = []
        self.attributes = None
        return
    def start_element(self, name, attributes):
        """ Handle opening element """
        self.attributes = attributes
        if name == "residue":
            obj = DefinitionResidue()
            self.curholder = obj
            self.curobj = obj
        elif name == "patch":
            obj = Patch()
            self.curholder = obj
            self.curobj = obj
        elif name == "atom":
            obj = DefinitionAtom()
            self.curatom = obj
            self.curobj = obj
        else:
            self.curelement = name
        return
    def end_element(self, name):
        """ Handle closing element """
        if name == "residue": # Complete Residue object
            residue = self.curholder
            if not isinstance(residue, DefinitionResidue):
                raise PDBInternalError("Internal error parsing XML!")
            resname = residue.name
            if resname == "":
                raise PDBInternalError("Residue name not set in XML!")
            else:
                self.map[resname] = residue
                self.curholder = None
                self.curobj = None
        elif name == "patch": # Complete patch object
            patch = self.curholder
            if not isinstance(patch, Patch):
                raise PDBInternalError("Internal error parsing XML!")
            patchname = patch.name
            if patchname == "":
                raise PDBInternalError("Residue name not set in XML!")
            else:
                self.patches.append(patch)
                self.curholder = None
                self.curobj = None
        elif name == "atom": # Complete atom object
            atom = self.curatom
            if not isinstance(atom, DefinitionAtom):
                raise PDBInternalError("Internal error parsing XML!")
            atomname = atom.name
            if atomname == "":
                raise PDBInternalError("Atom name not set in XML!")
            else:
                self.curholder.map[atomname] = atom
                self.curatom = None
                self.curobj = self.curholder
        else: # Just free the current element namespace
            self.curelement = ""
        return self.map
    def characters(self, text):
        """ Characters for the definition """
        if text.isspace():
            return
        # If this is a float, make it so
        try:
            value = float(str(text))
        except ValueError:
            value = str(text)
        # Special cases - lists and dictionaries
        if self.curelement == "bond":
            self.curobj.bonds.append(value)
        elif self.curelement == "dihedral":
            self.curobj.dihedrals.append(value)
        elif self.curelement == "altname":
            self.curholder.altnames[value] = self.curatom.name
        elif self.curelement == "remove":
            self.curobj.remove.append(value)
        else:
            setattr(self.curobj, self.curelement, value)
        return

class Definition:
    """ The Definition class contains the structured definitions found in the files and several
    mappings for easy access to the information. """
    def __init__(self):
        self.map = {}
        self.patches = {}
        handler = DefinitionHandler()
        sax.make_parser()
        for path in [AAPATH, NAPATH]:
            defpath = get_data_file(path)
            if defpath == "":
                raise PDBInternalError("%s not found!" % path)
            acid_file = open(defpath)
            sax.parseString(acid_file.read(), handler)
            acid_file.close()
            self.map.update(handler.map)
        # Now handle patches
        defpath = get_data_file(PATCHPATH)
        if defpath == "":
            raise PDBInternalError("%s not found!" % PATCHPATH)
        handler.map = {}
        patch_file = open(defpath)
        sax.parseString(patch_file.read(), handler)
        patch_file.close()
        # Apply specific patches to the reference object, allowing users to specify protonation
        # states in the PDB file
        for patch in handler.patches:
            if patch.newname != "":
                # Find all residues matching applyto
                resnames = list(self.map.keys())
                for name in resnames:
                    regexp = re.compile(patch.applyto).match(name)
                    if not regexp:
                        continue
                    newname = patch.newname.replace("*", name)
                    self.add_patch(patch, name, newname)
            # Either way, make sure the main patch name is available
            self.add_patch(patch, patch.applyto, patch.name)
    def add_patch(self, patch, refname, newname):
        """ Add a patch to a definition residue.
        Parameters
            patch: The patch object to add (Patch)
            refname   The name of the object to add the patch to (string)
            newname: The name of the new (patched) object (string) """
        try:
            aadef = self.map[refname] # The reference
            patch_residue = copy.deepcopy(aadef)
            # Add atoms from patch
            for atomname in patch.map:
                patch_residue.map[atomname] = patch.map[atomname]
                for bond in patch.map[atomname].bonds:
                    if bond not in patch_residue.map:
                        continue
                    if atomname not in patch_residue.map[bond].bonds:
                        patch_residue.map[bond].bonds.append(atomname)
            # Rename atoms as directed
            for key in patch.altnames:
                patch_residue.altnames[key] = patch.altnames[key]
            # Remove atoms as directed
            for remove in patch.remove:
                if not patch_residue.has_atom(remove):
                    continue
                removebonds = patch_residue.map[remove].bonds
                del patch_residue.map[remove]
                for bond in removebonds:
                    if remove in patch_residue.map[bond].bonds:
                        patch_residue.map[bond].bonds.remove(remove)
            # Add the new dihedrals
            for dihedral in patch.dihedrals:
                patch_residue.dihedrals.append(dihedral)
            # Point at the new reference
            self.map[newname] = patch_residue
            # Store the patch
            self.patches[newname] = patch
        except KeyError: # Just store the patch
            self.patches[newname] = patch

class Patch:
    """ Patch the definitionResidue class """
    def __init__(self):
        self.name = ""
        self.applyto = ""
        self.map = {}
        self.remove = []
        self.altnames = {}
        self.dihedrals = []
        self.newname = ""
    def __str__(self):
        text = "%s\n" % self.name
        text += "Apply to: %s\n" % self.applyto
        text += "Atoms to add: \n"
        for atom in self.map:
            text += "\t%s\n" % str(self.map[atom])
        text += "Atoms to remove: \n"
        for remove in self.remove:
            text += "\t%s\n" % remove
        text += "Alternate naming map: \n"
        text += "\t%s\n" % self.altnames
        return text

class DefinitionResidue(Residue):
    """ The DefinitionResidue class extends the Residue class to allow for a trimmed down
    initializing function. """
    def __init__(self):
        super().__init__()
        self.name = ""
        self.dihedrals = []
        self.map = {}
        self.altnames = {}
        self.dihedral_atoms = []
    def __str__(self):
        text = "%s\n" % self.name
        text += "Atoms: \n"
        for atom in self.map:
            text += "\t%s\n" % str(self.map[atom])
        text += "Dihedrals: \n"
        for dihedral in self.dihedrals:
            text += "\t%s\n" % dihedral
        text += "Alternate naming map: \n"
        text += "\t%s\n" % self.altnames
        return text
    def add_dihedral(self, atom):
        """ Add the atom to the list of dihedral bonds
        Parameters:
            atom: The atom to be added """
        self.dihedral_atoms.append(atom)
    def get_nearest_bonds(self, atomname):
        """ Get nearest bonds
        Parameters
            number: The number of bonds to get
        Returns
            bonds: A list of atomnames that are within three bonds of the atom and present in
            residue (list) """
        bonds = []
        lev2bonds = []
        atom = self.map[atomname]
        # Get directly bonded (length = 1) atoms
        for bondedatom in atom.bonds:
            if bondedatom not in bonds:
                bonds.append(bondedatom)
        # Get bonded atoms 2 bond lengths away
        for bondedatom in atom.bonds:
            for bond2 in self.map[bondedatom].bonds:
                if bond2 not in bonds and bond2 != atomname:
                    bonds.append(bond2)
                    lev2bonds.append(bond2)
        # Get bonded atoms 3 bond lengths away
        for lev2atom in lev2bonds:
            for bond3 in self.map[lev2atom].bonds:
                if bond3 not in bonds:
                    bonds.append(bond3)
        return bonds

class DefinitionAtom(Atom):
    """ A trimmed down version of the Atom class """
    def __init__(self, name=None, x=None, y=None, z=None):
        super().__init__()
        self.name = name
        self.x = x
        self.y = y
        self.z = z
        if name == None:
            self.name = ""
        if x == None:
            self.x = 0.0
        if y == None:
            self.y = 0.0
        if z == None:
            self.z = 0.0
        self.bonds = []
    def __str__(self):
        text = "%s: %.3f %.3f %.3f" % (self.name, self.x, self.y, self.z)
        for bond in self.bonds:
            text += " %s" % bond
        return text
    def isBackbone(self):
        """ Return true if atom name is in backbone, otherwise false """
        if self.name in BACKBONE:
            return 1
        else: return 0

