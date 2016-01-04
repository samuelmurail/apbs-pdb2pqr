""" This module contains the protein object used in PDB2PQR and associated methods

----------------------------

PDB2PQR -- An automated pipeline for the setup, execution, and analysis of Poisson-Boltzmann
electrostatics calculations

Copyright (c) 2002-2016, Jens Erik Nielsen, University College Dublin; Nathan A. Baker, Battelle
Memorial Institute, Developed at the Pacific Northwest National Laboratory, operated by Battelle
Memorial Institute, Pacific Northwest Division for the U.S. Department Energy.; Paul Czodrowski &
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
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

----------------------------
"""

__date__ = "2016-01-03"
__author__ = "Todd Dolinsky, Yong Huang, Nathan Baker"

import string

from .pdb import TER, ATOM, HETATM, END, MODEL
from .structures import Chain, Residue
from .nucleics import Nucleic

class Protein:
    """ The protein class represents the parsed PDB, and provides a hierarchy
    of information - each Protein contains a list of Chain objects as provided
    in the PDB file.  Each Chain then contains its associated list of Residue
    objects, and each Residue contains a list of Atom objects, completing the
    hierarchy.  """
    def __init__(self, pdblist, definition):
        """ Parameters
            pdblist: List of Classes of PDB lines as created """
        self.chainmap = {}
        self.chains = []
        self.residues = []
        self.referencemap = definition.map
        self.patchmap = definition.patches
        chain_dict = {}
        previous_atom = None
        residue = []
        num_models = 0
        num_chains = 1
        count = 0
        for record in pdblist: # Find number of chains
            if isinstance(record, TER):
                num_chains += 1
        for record in pdblist:
            if isinstance(record, ATOM) or isinstance(record, HETATM):
                if record.chain_id == "" and num_chains > 1 \
                and record.res_name not in ["WAT", "HOH"]:
                    # Assign a chain ID
                    record.chain_id = string.ascii_uppercase[count]
                chain_id = record.chain_id
                res_seq = record.res_seq
                insert_code = record.insert_code
                if previous_atom == None:
                    previous_atom = record
                if chain_id not in chain_dict:
                    my_chain = Chain(chain_id)
                    chain_dict[chain_id] = my_chain
                if res_seq != previous_atom.res_seq or \
                      insert_code != previous_atom.insert_code or \
                      chain_id != previous_atom.chain_id:
                    my_residue = self.create_residue(residue, previous_atom.res_name)
                    chain_dict[previous_atom.chain_id].add_residue(my_residue)
                    residue = []
                residue.append(record)
                previous_atom = record
            elif isinstance(record, END):
                my_residue = self.create_residue(residue, previous_atom.res_name)
                chain_dict[previous_atom.chain_id].add_residue(my_residue)
                residue = []
            elif isinstance(record, MODEL):
                num_models += 1
                if residue == []:
                    continue
                if num_models > 1:
                    my_residue = self.create_residue(residue, previous_atom.res_name)
                    chain_dict[previous_atom.chain_id].add_residue(my_residue)
                    break
            elif isinstance(record, TER):
                count += 1
        if residue != [] and num_models <= 1:
            my_residue = self.create_residue(residue, previous_atom.res_name)
            chain_dict[previous_atom.chain_id].add_residue(my_residue)
        # Keep a map for accessing chains via chain_id
        self.chainmap = chain_dict.copy()
        # Make a list for sequential ordering of chains
        if "" in chain_dict:
            chain_dict["ZZ"] = chain_dict[""]
            del chain_dict[""]
        keys = list(chain_dict.keys())
        keys.sort()
        for key in keys:
            self.chains.append(chain_dict[key])
        for chain in self.chains:
            for residue in chain.get_residues():
                self.residues.append(residue)
    def create_residue(self, residue, resname):
        """ Create a residue object.  If the resname is a known residue type, try to make that
        specific object, otherwise just make a standard residue object.

        Parameters
            residue:  A list of atoms (list)
            resname:  The name of the residue (string)

        Returns:
            residue:  The residue object (Residue) """
        try:
            refobj = self.referencemap[resname]
            if refobj.name != resname:
                # TODO - remove evals from this code!
                obj = "%s(residue, refobj)" % refobj.name
                residue = eval(obj)
                residue.reference = refobj
            else:
                # TODO - remove evals from this code!
                obj = "%s(residue, refobj)" % resname
                residue = eval(obj)
        except (KeyError, NameError):
            residue = Residue(residue)
        return residue
    def print_atoms(self, atomlist, chainflag=False, pdbfile=False):
        """ Get the text for the entire protein
        Parameters
            atomlist:  The list of atoms to include (list)
            chainflag: Flag whether to print chainid or not - Defaults to False
        Returns
            text:      The list of (stringed) atoms (list) """
        self.re_serialize()
        text = []
        currentchain_id = None
        for atom in atomlist:
            # Print the "TER" records between chains
            if currentchain_id == None:
                currentchain_id = atom.chain_id
            elif atom.chain_id != currentchain_id:
                currentchain_id = atom.chain_id
                text.append("TER\n")

            if pdbfile == True:
                text.append("%s\n" % atom.getPDBString())
            else:
                text.append("%s\n" % atom.getPQRString(chainflag=chainflag))
        text.append("TER\nEND")
        return text
    def re_serialize(self):
        """ Generate new serial numbers for atoms in the protein """
        count = 1
        for atom in self.get_atoms():
            atom.set("serial", count)
            count += 1
    def get_residues(self):
        """ Return the list of residues in the entire protein """
        return self.residues
    def num_residues(self):
        """ Get the number of residues for the entire protein (including multiple chains)
        Returns
            count:  Number of residues in the protein (int) """
        return len(self.get_residues())
    def num_atoms(self):
        """ Get the number of atoms for the entire protein(including multiple chains) """
        return len(self.get_atoms())
    def get_atoms(self):
        """ Return all Atom objects in list format.
        Returns
            atomlist:  List of Atom objects in the protein (list) """
        atomlist = []
        for chain in self.chains:
            for atom in chain.get_atoms():
                atomlist.append(atom)
        return atomlist
    def get_charge(self):
        """ Get the total charge on the protein
        NOTE: Since the misslist is used to identify incorrect charge assignments, this routine does
        not list the 3 and 5 termini of nucleic acid chains as having non-integer charge even though
        they are (correctly) non-integer.

        Returns:
            misslist: List of residues with non-integer charges (list)
            charge:   The total charge on the protein (float) """
        charge = 0.0
        misslist = []
        for chain in self.chains:
            for residue in chain.get("residues"):
                rescharge = residue.get_charge()
                charge += rescharge
                if isinstance(residue, Nucleic):
                    if residue.is3term or residue.is5term:
                        continue
                if float("%i" % rescharge) != rescharge:
                    misslist.append(residue)
        return misslist, charge
    def get_chains(self):
        """ Get the chains object
        Returns
            chains: The list of chains in the protein (chain) """
        return self.chains
    def get_summary(self):
        """ Get a summary of the protein object """
        output = []
        for chain in self.chains:
            output.append(chain.get_summary())
        return ' '.join(output)
