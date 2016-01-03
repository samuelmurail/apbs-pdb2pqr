""" Structures for PDB2PQR

    This module contains the structure objects used in PDB2PQR and their associated methods.

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
    * Neither the names of University College Dublin, Battelle Memorial Institute,Pacific Northwest
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

BACKBONE = ["N", "CA", "C", "O", "O2", "HA", "HN", "H", "tN"]

from . import pdb
from . import utilities
from . import quatfit
from .errors import PDBInternalError

class Chain:
    """ The chain class contains information about each chain within a given Protein object. """
    def __init__(self, chain_id):
        """ chain_id: The chain_id for this chain as denoted in the PDB file (string) """
        self.chain_id = chain_id
        self.residues = []
    def get_member(self, name):
        """ Get a member of the Chain class
        Parameters
        * name: The name of the member; possible values:
            * ID: The ID of the chain
            * Residues: The list of residues within the Chain
        Returns
        * item: The value of the member """
        if name == "atoms":
            self.get_atoms()
        else:
            try:
                item = getattr(self, name)
                return item
            except AttributeError:
                message = 'Unable to get object "%s" in class Chain' % name
                raise PDBInternalError(message)
    def add_residue(self, residue):
        """ Add a residue to the chain
        Parameters
        * residue: The residue to be added (Residue)
        """
        self.residues.append(residue)
    def num_residues(self):
        """ Get the number of residues for the chain. """
        return len(self.residues)
    def renumber_residues(self):
        """ Renumber Atoms based on actual Residue number and not PDB numbering. """
        count = 1
        for residue in self.residues:
            residue.setResSeq(count)
            count += 1
    def num_atoms(self):
        """ Get the number of atoms for the chain. """
        return len(self.get_atoms())
    def get_residues(self):
        """ Return a list of Residue objects in this chain"""
        return self.residues
    def get_atoms(self):
        """ Return a list of Atom objects contained in this chain """
        atomlist = []
        for residue in self.residues:
            res_atom_list = residue.get("atoms")
            for atom in res_atom_list:
                atomlist.append(atom)
        return atomlist
    def get_summary(self):
        """ Return a summary string for the object """
        # TODO - why isn't this just the __str__ member?
        output = []
        for residue in self.residues:
            output.append(residue.letter_code())
        return ''.join(output)

class Residue:
    """ The residue class contains a list of Atom objects associated with that residue and other
    helper functions. """
    def __init__(self, atoms):
        """ atoms: A list of Atom objects to be stored in this class (list) """
        sample_atom = atoms[-1]
        self.atoms = []
        self.name = sample_atom.res_name
        self.chain_id = sample_atom.chain_id
        self.res_seq = sample_atom.res_seq
        self.insert_code = sample_atom.insert_code
        self.map = {}
        self.naname = None
        atomclass = ""
        for atom in atoms:
            if isinstance(atom, pdb.ATOM):
                atomclass = "ATOM"
            elif isinstance(atom, pdb.HETATM):
                atomclass = "HETATM"
            atom = Atom(atom, atomclass, self)
            atomname = atom.get("name")
            if atomname not in self.map:
                self.add_atom(atom)
            else: # Don't add duplicate atom
                oldatom = self.get_atom(atomname)
                oldatom.set("alt_loc", "")
        if self.name == "HOH":
            self.name = "WAT"
            for atom in self.atoms:
                atom.set("res_name", "WAT")
    def __str__(self):
        text = "%s %s %i%s" % (self.name, self.chain_id, self.res_seq, self.insert_code)
        return text
    def update_terminus_status(self):
        """ Update the is_n_terminus and is_c_terminus flags"""
        if self.is_n_terminus:
            count = 0
            atoms = ['H', 'H2', 'H3']
            for atom in atoms:
                for atom2 in self.atoms:
                    atomname = atom2.get('name')
                    if atom == atomname:
                        count = count+1
            self.is_n_terminus = count
        if self.is_c_terminus:
            self.is_c_terminus = None
            for atom in self.atoms:
                atomname = atom.get('name')
                if atomname == 'HO':
                    # TODO - this is insanity... why is a Boolean being updated as an int?
                    self.is_c_terminus = 2
                    break
            if not self.is_c_terminus:
                self.is_c_terminus = 1
        return
    def num_atoms(self):
        """ Get the number of atoms for the residue"""
        return len(self.atoms)
    def setResSeq(self, value):
        """ Set the atom field res_seq to a certain value andchange the residue's information.  The
        icode field is no longer useful.
        Parameters
            value:  The new value of res_seq (int) """
        self.insert_code = ""
        self.res_seq = value
        for atom in self.atoms:
            atom.set("res_seq", value)
    def set_chain_id(self, value):
        """ Set the chain_id field to a certain value """
        self.chain_id = value
        for atom in self.atoms:
            atom.set("chain_id", value)
    def add_atom(self, atom):
        """ Add the atom object to the residue.
        Parameters:
            atom - The object to be added (ATOM)
        """
        self.atoms.append(atom)
        self.map[atom.get("name")] = atom
    def remove_atom(self, atomname):
        """ Remove an atom from the residue object.
        Parameters
            atomname: The name of the atom to be removed (string)
        """
        atom = self.map[atomname]
        bonds = atom.bonds
        del self.map[atomname]
        self.atoms.remove(atom)
        for bondatom in bonds:
            if atom in bondatom.bonds:
                bondatom.bonds.remove(atom)
        del atom
    def rename_atom(self, oldname, newname):
        """ Rename an atom to a new name
        Parameters
            oldname: The old atom name (string)
            newname: The new atom name (string)
        """
        atom = self.map[oldname]
        atom.set("name", newname)
        self.map[newname] = atom
        del self.map[oldname]
    def create_atom(self, name, newcoords, atom_type):
        """ Add a new atom object to the residue. Uses an atom currently in the residue to seed the
        new atom object, then replaces the coordinates and name accordingly.
        Parameters
            name: The name of the new atom (string)
            newcoords: The x,y,z coordinates of the new atom (list)
            type: The type of atom, ATOM or HETATM """
        oldatom = self.atoms[0]
        newatom = Atom(oldatom, atom_type, self)
        newatom.set("x", newcoords[0])
        newatom.set("y", newcoords[1])
        newatom.set("z", newcoords[2])
        newatom.set("name", name)
        newatom.set("occupancy", 1.00)
        newatom.set("temperature_factor", 0.00)
        self.add_atom(newatom)
    def addMissing(self, value):
        """ Add the value to the list of missing atoms
        Parameters
            value: The name of the missing atom (string)"""
        self.missing.append(value)
    def get_atom(self, name):
        """ Retrieve an atom from the mapping
        Parameters
            resname: The name of the residue to retrieve (string)"""
        return self.map.get(name)
    def get_atoms(self):
        """ Get the atoms """
        return self.atoms
    def has_atom(self, name):
        """ Does it have this atom? """
        return name in self.map
    def getCharge(self):
        """ Get the total charge of the residue.  In order to get rid of floating point rounding
        error, do the string transformation.
        """
        charge = (atom.ffcharge for atom in self.atoms if atom.ffcharge)
        charge = sum(charge)
        charge = float("%.4f" % charge)
        return charge
    def renameResidue(self, name):
        """ Rename a given residue
        Parameters
            name:       The new name of the residue """
        self.name = name
        for atom in self.atoms:
            atom.res_name = name
    def rotateTetrahedral(self, atom1, atom2, angle):
        """ Rotate about the atom1-atom2 bond by a given angle All atoms connected to atom2 will
        rotate.
        Parameters:
            atom1:  The first atom of the bond to rotate about (atom)
            atom2:  The second atom of the bond to rotate about (atom)
            angle:  The number of degrees to rotate (float) """
        moveatoms = []
        movecoords = []
        initcoords = subtract(atom2.getCoords(), atom1.getCoords())
        for atom in atom2.bonds:
            if atom == atom1: continue
            moveatoms.append(atom)
            movecoords.append(subtract(atom.getCoords(), atom1.getCoords()))
        newcoords = qchichange(initcoords, movecoords, angle)
        for i in range(len(moveatoms)):
            atom = moveatoms[i]
            x = (newcoords[i][0] + atom1.get("x"))
            y = (newcoords[i][1] + atom1.get("y"))
            z = (newcoords[i][2] + atom1.get("z"))
            atom.set("x", x)
            atom.set("y", y)
            atom.set("z", z)
    def set_donors_and_acceptors(self):
        """ Set the donors and acceptors within the residue """
        if not hasattr(self, "reference"): return
        for atom in self.get_atoms():
            atomname = atom.get("name")
            resname = self.name
            # TODO - are 0s and 1s being used for True and False again?
            atom.set("hdonor", 0)
            atom.set("hacceptor", 0)
            if atomname.startswith("N"):
                bonded = 0
                for bondedatom in atom.bonds:
                    if bondedatom.isHydrogen():
                        atom.set("hdonor", 1)
                        bonded = 1
                        break
                if not bonded and self.reference.name == "HIS":
                    atom.set("hacceptor", 1)
            elif atomname.startswith("O") or \
                 (atomname.startswith("S") and self.reference.name == "CYS"):
                atom.set("hacceptor", 1)
                for bondedatom in atom.bonds:
                    if bondedatom.isHydrogen():
                        atom.set("hdonor", 1)
                        break
    def reorder(self):
        """ Reorder the atoms to start with N, CA, C, O if they exist """
        templist = []
        if self.has_atom("N"): templist.append(self.get_atom("N"))
        if self.has_atom("CA"): templist.append(self.get_atom("CA"))
        if self.has_atom("C"): templist.append(self.get_atom("C"))
        if self.has_atom("O"): templist.append(self.get_atom("O"))
        for atom in self.atoms:
            if atom.name not in ["N", "CA", "C", "O"]:
                templist.append(atom)
        self.atoms = templist[:]
    def letter_code(self):
        return 'X'

class Atom(pdb.ATOM):
    """ The Atom class inherits off the ATOM object in pdb.py.  It is used for adding fields not
    found in the pdb that may be useful for analysis. Also simplifies code by combining ATOM and
    HETATM objects into a single class. """
    def __init__(self, atom, type, residue):
        """ Parameters
        atom: The original ATOM object (ATOM)
        type: Either ATOM or HETATM (string)
        residue: A pointer back to the parent residue object (Residue)
        """
        if type == "ATOM" or type == "HETATM":
            self.type = type
        else:
            raise PDBInternalError("Invalid atom type %s (Atom Class IN structures.py)!" % type)
        self.serial = atom.serial
        self.name = atom.name
        self.alt_loc = atom.alt_loc
        self.res_name = atom.res_name
        self.chain_id = atom.chain_id
        self.res_seq = atom.res_seq
        self.insert_code = atom.insert_code
        self.x = atom.x
        self.y = atom.y
        self.z = atom.z
        self.occupancy = atom.occupancy
        self.temperature_factor = atom.temperature_factor
        self.segment_id = atom.segment_id
        self.element = atom.element
        self.charge = atom.charge
        self.bonds = []
        self.reference = None
        self.residue = residue
        self.radius = None
        self.ffcharge = None
        self.hdonor = 0
        self.hacceptor = 0
        self.cell = None
        self.added = 0
        self.optimizeable = 0
        self.refdistance = 0
        self.id = None
        self.mol2charge = None
        if hasattr(atom, 'mol2charge'):
            self.mol2charge = atom.mol2charge
    def getCommonStringRep(self, chainflag=False):
        """ Returns a string of the common column of the new atom type. Uses the ATOM string output
        but changes the first field to either by ATOM or HETATM as necessary. This is used to create
        the output for pqr and pdb files!
        Returns string with ATOM/HETATM field set appropriately """
        outstr = ""
        tstr = self.type
        outstr += string.ljust(tstr, 6)[:6]
        tstr = "%d" % self.serial
        outstr += string.rjust(tstr, 5)[:5]
        outstr += " "
        tstr = self.name
        if len(tstr) == 4 or len(tstr.strip("FLIP")) == 4:
            outstr += string.ljust(tstr, 4)[:4]
        else:
            outstr += " " + string.ljust(tstr, 3)[:3]
        tstr = self.res_name
        if len(tstr) == 4:
            outstr += string.ljust(tstr, 4)[:4]
        else:
            outstr += " " + string.ljust(tstr, 3)[:3]
        outstr += " "
        if chainflag:
            tstr = self.chain_id
        else:
            tstr = ''
        outstr += string.ljust(tstr, 1)[:1]
        tstr = "%d" % self.res_seq
        outstr += string.rjust(tstr, 4)[:4]
        if self.insert_code != "":
            outstr += "%s   " % self.insert_code
        else:
            outstr += "    "
        tstr = "%8.3f" % self.x
        outstr += string.ljust(tstr, 8)[:8]
        tstr = "%8.3f" % self.y
        outstr += string.ljust(tstr, 8)[:8]
        tstr = "%8.3f" % self.z
        outstr += string.ljust(tstr, 8)[:8]
        return outstr
    def __str__(self):
        """ Returns a string of the new atom type.  Uses the ATOM string output but changes the
        first field to either by ATOM or HETATM as necessary. This is used to create the output for
        pqr files! Returns string with ATOM/HETATM field set appropriately """
        return self.getPQRString()
    def getPQRString(self, chainflag=False):
        """ Returns a string of the new atom type.  Uses the ATOM string output but changes the
        first field to either by ATOM or HETATM as necessary. This is used to create the output for
        pqr files! Returns string with ATOM/HETATM field set appropriately """
        outstr = self.getCommonStringRep(chainflag=chainflag)
        if self.ffcharge != None:
            ffcharge = "%.4f" % self.ffcharge
        else:
            ffcharge = "0.0000"
        outstr += string.rjust(ffcharge, 8)[:8]
        if self.radius != None:
            ffradius = "%.4f" % self.radius
        else:
            ffradius = "0.0000"
        outstr += string.rjust(ffradius, 7)[:7]
        return outstr
    def getPDBString(self):
        """ Returns a string of the new atom type.  Uses the ATOM string output but changes the
        first field to either by ATOM or HETATM as necessary. This is for the pdb representation of
        the atom. The propka30 module depends on this being correct. Returns string with ATOM/HETATM
        field set appropriately """
        outstr = self.getCommonStringRep(chainflag=True)
        tstr = "%6.2f" % self.occupancy
        outstr += string.ljust(tstr, 6)[:6]
        tstr = "%6.2f" % self.temperature_factor
        outstr += string.rjust(tstr, 6)[:6]
        #padding between temp factor and segment_id
        outstr += ' ' * 7
        tstr = self.segment_id
        outstr += string.ljust(tstr, 4)[:4]
        tstr = self.element
        outstr += string.ljust(tstr, 2)[:2]
        tstr = str(self.charge)
        outstr += string.ljust(tstr, 2)[:2]
        return outstr
    def getCoords(self):
        """ Return the x,y,z coordinates of the atom in list form """
        return [self.x, self.y, self.z]
    def addBond(self, bondedatom):
        """ Add a bond to the list of bonds
        Parameters:
            bondedatom: The atom to bond to (Atom) """
        self.bonds.append(bondedatom)
    def isHydrogen(self):
        """ Is this atom a hydrogen? """
        return self.name[0] == "H"
    def isBackbone(self):
        """ Is this atom in the backbone? """
        return self.name in BACKBONE
    def hasReference(self):
        """ Determine if the atom object has a reference object or not. All known atoms should have
        reference objects. """
        return self.reference != None
