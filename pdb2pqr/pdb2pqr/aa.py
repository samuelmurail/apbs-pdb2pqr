""" Amino Acid Structures for PDB2PQR

    This module contains the base amino acid structures for pdb2pqr.

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

from structures import Atom, Residue
from .errors import PDBInputError

class Amino(Residue):
    """ This class provides standard features of the amino acids listed below
    Parameters
    * atoms:  A list of Atom objects to be stored in this class (list)
    * ref:  The reference object for the amino acid.  Used to convert from the alternate naming
      scheme to the main naming scheme. """
    def __init__(self, atoms, ref):
        sample_atom = atoms[-1]
        self.atoms = []
        self.name = sample_atom.res_name
        self.chain_id = sample_atom.chain_id
        self.res_seq = sample_atom.res_seq
        self.insert_code = sample_atom.insert_code
        self.ffname = self.name
        self.map = {}
        self.dihedrals = []
        self.patches = []
        self.peptide_c = None
        self.peptide_n = None
        self.is_n_terminus = False
        self.is_c_terminus = False
        self.is5term = False
        self.is3term = False
        self.missing = []
        self.reference = ref
        self.fixed = 0
        self.stateboolean = {}
        for atom in atoms:
            if atom.name in ref.altnames: # Rename atoms
                atom.name = ref.altnames[atom.name]
            if atom.name not in self.map:
                atom_ = Atom(atom, "ATOM", self)
                self.add_atom(atom_)
    def create_atom(self, atomname, newcoords):
        """ Create an atom.  Override the generic residue's version of create_atom().
        Parameters
        * atomname:  The name of the atom (string)
        * newcoords: The coordinates of the atom (list). """
        oldatom = self.atoms[0]
        newatom = Atom(oldatom, "ATOM", self)
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
    def add_dihedral_angle(self, value):
        """ Add the value to the list of chi angles
        Parameters
        * value: The value to be added (float)"""
        self.dihedrals.append(value)
    def set_state(self):
        """ Set the name to use for the forcefield based on the current state.  Uses N* and C* for
        termini. """
        if self.is_n_terminus:
            if "NEUTRAL-NTERM" in self.patches:
                self.ffname = "NEUTRAL-N%s" % self.ffname
            else:
                self.ffname = "N%s" % self.ffname
        elif self.is_c_terminus:
            if "NEUTRAL-CTERM" in self.patches:
                self.ffname = "NEUTRAL-C%s" % self.ffname
            else:
                self.ffname = "C%s" % self.ffname
    def letter_code(self):
        """ Return one-letter amino acid abbreviation """
        raise NotImplementedError

class ALA(Amino):
    """ This class gives data about the Alanine object, and inherits off the base residue class. """
    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref
    def letter_code(self):
        return 'A'

class ARG(Amino):
    """ This class gives data about the Arginine object, and inherits off the base residue
    class. """
    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref
    def letter_code(self):
        return 'R'
    def set_state(self):
        """ Set the name to use for the forcefield based on the current state. """
        if "AR0" in self.patches or self.name == "AR0":
            self.ffname = "AR0"
        Amino.set_state(self)

class ASN(Amino):
    """ This class gives data about the Asparagine object, and inherits off the base residue
    class. """
    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref
    def letter_code(self):
        return 'N'

class ASP(Amino):
    """ This class gives data about the Aspartic Acid object, and inherits off the base residue
    class. """
    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref
    def letter_code(self):
        return 'D'
    def set_state(self):
        """ Set the name to use for the forcefield based on the current state. """
        if "ASH" in self.patches or self.name == "ASH":
            self.ffname = "ASH"
        Amino.set_state(self)

class CYS(Amino):
    """ This class gives data about the Cysteine object, and inherits off the base residue
    class. """
    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref
        self.disulfide_bonded = 0
        self.disulfide_bonded_partner = None
    def letter_code(self):
        return 'C'
    def set_state(self):
        """ Set the state of the CYS object.  If SS-bonded, use CYX. If negatively charged, use CYM.
        If HG is not present, use CYX. """
        if "CYX" in self.patches or self.name == "CYX":
            self.ffname = "CYX"
        elif self.disulfide_bonded:
            self.ffname = "CYX"
        elif "CYM" in self.patches or self.name == "CYM":
            self.ffname = "CYM"
        elif not self.has_atom("HG"):
            self.ffname = "CYX"
        Amino.set_state(self)

class GLN(Amino):
    """ This class gives data about the Glutamine object, and inherits off the base residue
    class. """
    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref
    def letter_code(self):
        return 'Q'

class GLU(Amino):
    """ This class gives data about the Glutamic Acid object, and inherits off the base residue
    class. """
    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref
    def letter_code(self):
        return 'E'
    def set_state(self):
        """  Set the name to use for the forcefield based on the current state."""
        if "GLH" in self.patches or self.name == "GLH":
            self.ffname = "GLH"
        Amino.set_state(self)

class GLY(Amino):
    """ This class gives data about the Glycine object, and inherits off the base residue class. """
    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref
    def letter_code(self):
        return 'G'

class HIS(Amino):
    """ This class gives data about the Histidine object, and inherits off the base residue
    class. """
    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref
    def letter_code(self):
        return 'H'
    def set_state(self):
        """ Histidines are a special case due to the presence of several different forms.  This
        function sets all non-positive incarnations of HIS to neutral HIS by checking to see if
        optimization removed hacceptor or hdonor flags.  Otherwise HID is used as the default. """
        if "HIP" not in self.patches and self.name not in ["HIP", "HSP"]:
            if self.get_atom("ND1").hdonor and not \
                   self.get_atom("ND1").hacceptor:
                if self.has_atom("HE2"):
                    self.remove_atom("HE2")
            elif self.get_atom("NE2").hdonor and not \
                     self.get_atom("NE2").hacceptor:
                if self.has_atom("HD1"):
                    self.remove_atom("HD1")
            elif self.get_atom("ND1").hacceptor and not \
                     self.get_atom("ND1").hdonor:
                if self.has_atom("HD1"):
                    self.remove_atom("HD1")
            else: # Default to HID
                if self.has_atom("HE2"):
                    self.remove_atom("HE2")
        if self.has_atom("HD1") and self.has_atom("HE2"):
            self.ffname = "HIP"
        elif self.has_atom("HD1"):
            self.ffname = "HID"
        elif self.has_atom("HE2"):
            self.ffname = "HIE"
        else:
            raise PDBInputError("Invalid type for %s! Missing both HD1 and HE2 atoms."
                                " If you receive this error while using the --assign-only"
                                " option you can only resolve it by adding HD1, HE2 or both to"
                                " this residue." % str(self))
        Amino.set_state(self)

class ILE(Amino):
    """ This class gives data about the Isoleucine object, and inherits off the base residue
    class. """
    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref
    def letter_code(self):
        return 'I'

class LEU(Amino):
    """ This class gives data about the Leucine object, and inherits off the base residue class. """
    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref
    def letter_code(self):
        return 'L'

class LYS(Amino):
    """ This class gives data about the Lysine object, and inherits off the base residue class. """
    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref
    def letter_code(self):
        return 'K'
    def set_state(self):
        """ Determine if this is LYN or not """
        if "LYN" in self.patches or self.name == "LYN":
            self.ffname = "LYN"
        Amino.set_state(self)

class MET(Amino):
    """ This class gives data about the Methionine object, and inherits off the base residue
    class. """
    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref
    def letter_code(self):
        return 'M'

class PHE(Amino):
    """ This class gives data about the Phenylalanine object, and inherits off the base residue
    class. """
    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref
    def letter_code(self):
        return 'F'

class PRO(Amino):
    """ This class gives data about the Proline object, and inherits off the base residue class. """
    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref
    def letter_code(self):
        return 'P'
    def set_state(self):
        """ Set the name to use for the forcefield based on the current state. Uses N* and C* for
        termini. """
        if self.is_n_terminus:
            self.ffname = "N%s" % self.ffname
        elif self.is_c_terminus:
            if "NEUTRAL-CTERM" in self.patches:
                self.ffname = "NEUTRAL-C%s" % self.ffname
            else:
                self.ffname = "C%s" % self.ffname

class SER(Amino):
    """ This class gives data about the Serine object, and inherits off the base residue class. """
    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref
    def letter_code(self):
        return 'S'

class THR(Amino):
    """ This class gives data about the Threonine object, and inherits off the base residue
    class. """
    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref
    def letter_code(self):
        return 'T'

class TRP(Amino):
    """ This class gives data about the Tryptophan object, and inherits off the base residue
    class. """
    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref
    def letter_code(self):
        return 'W'

class TYR(Amino):
    """ This class gives data about the Tyrosine object, and inherits off the base residue
    class. """
    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref
    def letter_code(self):
        return 'Y'
    def set_state(self):
        """ See if the TYR is negative or not """
        if "TYM" in self.patches or self.name == "TYM":
            self.ffname = "TYM"
        Amino.set_state(self)

class VAL(Amino):
    """ This class gives data about the Valine object, and inherits off the base residue class. """
    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref
    def letter_code(self):
        return 'V'

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
