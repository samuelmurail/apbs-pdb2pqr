""" Routines for PDB2PQR

    This module contains the protein object used in PDB2PQR and methods used to correct, analyze,
    and optimize that protein.

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

# TODO - pylint hates this module... it desperately needs cleaning up.

__date__ = "2016-01-02"
__author__ = "Jens Erik Nielsen, Todd Dolinsky, Yong Huang"

CELL_SIZE = 2
BUMP_DIST = 2.0
BUMP_HDIST = 1.5
BUMP_HYDROGEN_SIZE = 0.5
BUMP_HEAVY_SIZE = 1.0
BONDED_SS_LIMIT = 2.5
PEPTIDE_DIST = 1.7
REPAIR_LIMIT = 10
ANGLE_STEPS = 72
ANGLE_STEP_SIZE = float(360 // ANGLE_STEPS)
ANGLE_TEST_COUNT = 10
EPSILON = 0.0000001

AMINO_ACIDS = ["ALA", "ARG", "ASH", "ASN", "ASP", "CYS", "CYM", "GLN", "GLU", "GLH", "GLY",
               "HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "ILE", "LEU", "LYS", "LYN",
               "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "TYM", "VAL"]
NUCLEIC_ACIDS = ["A", "A5", "A3", "C", "C5", "C3", "G", "G5", "G3", "T", "T5", "T3", "U",
                 "U5", "U3", "RA", "RG", "RC", "RU", "DA", "DG", "DC", "DT"]

import copy
import sys
import string
from io import StringIO
from pprint import pformat

from .errors import PDBInputError, PDBInternalError, PDB2PKAError
from . import molecules as mols
from . import aminos as aa
from . import nucleics as na
from .structures import Chain
from .utilities import distance, get_dihedral, shortest_path, subtract
from .quatfit import find_coordinates, qchichange

class Routines:
    """ This is a container class for storing variables and running different functionality in
    PDB2PKA/PROPKA/etc.  """
    def __init__(self, protein, verbose, definition=None):
        """ Parameters
        protein:  The protein to run PDB2PQR on (Protein)
        verbose:  A flag to determine whether to write to stdout """
        self.protein = protein
        self.definition = definition
        self.aadef = None
        self.verbose = verbose
        self.warnings = []
        self.cells = {}
        if definition != None:
            self.aadef = definition.getAA()
            self.nadef = definition.getNA()
    def write(self, message, indent=0):
        """ Write a message to stdout for debugging if verbose
        Parameters
            message: The message to write (string)
            indent : The indent level (int, default=0) """
        out = ""
        if self.verbose:
            for _ in range(indent):
                out += "\t"
            out += message
            sys.stdout.write(out)
    def get_warnings(self):
        """ Get all warnings generated from routines """
        return self.warnings
    def apply_name_scheme(self, forcefield):
        """ Apply the naming scheme of the given forcefield to the atoms within the protein
        Parameters
            forcefield: The forcefield object (forcefield) """
        self.write("Applying the naming scheme to the protein...")
        for residue in self.protein.get_residues():
            if isinstance(residue, (aa.Amino, mols.WAT, na.Nucleic)):
                resname = residue.ffname
            else:
                resname = residue.name
            for atom in residue.get_atoms():
                rname, aname = forcefield.get_names(resname, atom.name)
                if resname not in ['LIG', 'WAT', 'ACE', 'NME'] and rname != None:
                    try:
                        if (residue.is_n_terminus or residue.is_c_terminus) \
                        and rname != residue.name:
                            rname = residue.name
                    except AttributeError:
                        pass
                if aname != None and rname != None:
                    atom.res_name = rname
                    atom.name = aname
        self.write("Done.\n")
    def apply_forcefield(self, forcefield):
        """ Apply the forcefield to the atoms within the protein
        Parameters
            forcefield: The forcefield object (forcefield)
        Returns
            hitlist:    A list of atoms that were found in the forcefield (list)
            misslist:   A list of atoms that were not found in the forcefield (list) """
        self.write("Applying the forcefield to the protein...")
        misslist = []
        hitlist = []
        for residue in self.protein.get_residues():
            if isinstance(residue, (aa.Amino, mols.WAT, na.Nucleic)):
                resname = residue.ffname
            else:
                resname = residue.name
            for atom in residue.get_atoms():
                atomname = atom.get("name")
                charge, radius = forcefield.get_params(resname, atomname)
                if charge != None and radius != None:
                    atom.set("ffcharge", charge)
                    atom.set("radius", radius)
                    hitlist.append(atom)
                else:
                    misslist.append(atom)
        self.write("Done.\n")
        return hitlist, misslist
    def update_residue_types(self):
        """ Find the type of residue as notated in the Amino Acid definition """
        self.write("Updating Residue Types... ")
        for chain in self.protein.getChains():
            for residue in chain.get("residues"):
                name = residue.get("name")
                if name in AMINO_ACIDS:
                    residue.set("type", 1)
                elif name == "WAT":
                    residue.set("type", 3)
                elif name in NUCLEIC_ACIDS:
                    residue.set("type", 4)
                else: # Residue is a ligand or unknown
                    residue.set("type", 2)
        self.write("Done\n")
    def update_disulfide_bridges(self):
        """ Check for SS-bridge partners, and if present, set appropriate partners """
        self.write("Updating SS bridges...\n")
        sg_partners = {}
        for residue in self.protein.get_residues():
            if isinstance(residue, aa.CYS):
                atom = residue.get_atom("SG")
                if atom != None:
                    sg_partners[atom] = []
        for atom in sg_partners:
            for partner in sg_partners:
                if atom == partner or sg_partners[atom] != []:
                    continue
                dist = distance(atom.getCoords(), partner.getCoords())
                if dist < BONDED_SS_LIMIT:
                    sg_partners[atom].append(partner)
                    sg_partners[partner].append(atom)
        for atom in sg_partners:
            res1 = atom.get("residue")
            numpartners = len(sg_partners[atom])
            if numpartners == 1:
                partner = sg_partners[atom][0]
                res2 = partner.get("residue")
                res1.set("disulfide_bonded", 1)
                res1.set("disulfide_bonded_partner", partner)
                self.apply_patch("CYX", res1)
                self.write("%s - %s\n" % (res1, res2), 1)
            elif numpartners > 1:
                error = "WARNING: %s has multiple potential " % res1
                error += "SS-bridge partners\n"
                self.write(error, 1)
                self.warnings.append(error)
            elif numpartners == 0:
                self.write("%s is a free cysteine\n" % res1, 1)
        self.write("Done.\n")
    def update_internal_linked_bonds(self):
        """ Update the internal bonding network using the reference objects in each atom. """
        for residue in self.protein.get_residues():
            if isinstance(residue, (aa.Amino, mols.WAT, na.Nucleic)):
                for atom in residue.get_atoms():
                    if not atom.hasReference():
                        continue
                    for bond in atom.reference.bonds:
                        if not residue.has_atom(bond):
                            continue
                        bondatom = residue.get_atom(bond)
                        if bondatom not in atom.bonds:
                            atom.addBond(bondatom)
    def update_bonds(self):
        """ Update the bonding network of the protein.  This happens in 3 steps:
        1.  Applying the PEPTIDE patch to all Amino residues so as to add reference for the N(i+1)
            and C(i-1) atoms
        2.  UpdateInternalinked_bonds for inter-residue linking
        3.  Set the links to the N(i+1) and C(i-1) atoms """
        for residue in self.protein.get_residues():
            if isinstance(residue, aa.Amino):
                if residue.is_n_terminus or residue.is_c_terminus:
                    continue
                else:
                    self.apply_patch("PEPTIDE", residue)
        self.update_internal_linked_bonds()
        for chain in self.protein.getChains():
            for i in range(chain.num_residues() - 1):
                res1 = chain.residues[i]
                res2 = chain.residues[i + 1]
                if not isinstance(res1, aa.Amino) or not isinstance(res2, aa.Amino):
                    continue
                atom1 = res1.get_atom("C")
                atom2 = res2.get_atom("N")
                if atom1 != None:
                    res2.peptide_c = atom1
                if atom2 != None:
                    res1.peptide_n = atom2
                if atom1 == None or atom2 == None:
                    continue
                if distance(atom1.getCoords(), atom2.getCoords()) > PEPTIDE_DIST:
                    text = "Gap in backbone detected between %s and %s!\n" % \
                           (res1, res2)
                    self.write(text, 1)
                    self.warnings.append(text)
                    res2.peptide_c = None
                    res1.peptide_n = None
    def apply_patch(self, patchname, residue):
        """ Apply a patch to the given residue.  This is one of the key functions in PDB2PQR.  A
        similar function appears in definitions.py - that version is needed for residue level
        substitutions so certain protonation states (i.e. CYM, HSE) are detectable on input.

        This version looks up the particular patch name in the patchmap stored in the protein, and
        then applies the various commands to the reference and actual residue structures.

        Parameters
            patchname:  The name of the patch (string)
            residue:    The residue to apply the patch to (residue) """
        if patchname not in self.protein.patchmap:
            raise PDBInternalError("Unable to find patch %s!" % patchname)

        self.write('PATCH INFO: %s patched with %s\n' % (residue, patchname), 1)
        # Make a copy of the reference, i.e. a new reference for
        # this patch.  Two examples:
        #     PEPTIDE is a special case, as it applies to
        #             every residue.
        #     CTERM only applies to one specific residue, so a
        #             deep copy is used.
        if patchname == "PEPTIDE":
            newreference = residue.reference
        else:
            newreference = copy.deepcopy(residue.reference)
        patch = self.protein.patchmap[patchname]
        # Add atoms from patch
        for atomname in patch.map:
            newreference.map[atomname] = patch.map[atomname]
            for bond in patch.map[atomname].bonds:
                if bond not in newreference.map:
                    continue
                if atomname not in newreference.map[bond].bonds:
                    newreference.map[bond].bonds.append(atomname)
        # Remove atoms as directed by patch
        for remove in patch.remove:
            if remove in residue.map:
                residue.remove_atom(remove)
            if remove not in newreference.map:
                continue
            removebonds = newreference.map[remove].bonds
            del newreference.map[remove]
            for bond in removebonds:
                index = newreference.map[bond].bonds.index(remove)
                del newreference.map[bond].bonds[index]
        # Add the new dihedrals
        for dihedral in patch.dihedrals:
            newreference.dihedrals.append(dihedral)
        # Point at the new reference
        residue.reference = newreference
        residue.patches.append(patchname)
        # Rename atoms as directed by patch
        for atom in residue.get_atoms():
            if atom.name in patch.altnames:
                residue.rename_atom(atom.name, patch.altnames[atom.name])
        # Replace each atom's reference with the new one
        for atomname in residue.map:
            if newreference.has_atom(atomname):
                atom = residue.get_atom(atomname)
                atom.reference = newreference.map[atomname]
    def set_states(self):
        """ Set the state of each residue.  This is the last step before assigning the forcefield,
        but is necessary so as to distinguish between various protonation states.
        See aa.py for residue-specific functions. """
        for residue in self.protein.get_residues():
            if isinstance(residue, (aa.Amino, na.Nucleic)):
                residue.set_state()
    def assign_termini(self, chain, neutraln=False, neutralc=False):
        """ Assign the termini for the given chain by looking at the start and end residues. """
        if len(chain.residues) == 0:
            text = "Error: chain \"%s\" has 0 residues!" % chain.chain_id
            raise PDBInputError(text)
        # Set the N-Terminus/ 5' Terminus
        res0 = chain.residues[0]
        if isinstance(res0, aa.Amino):
            res0.set("is_n_terminus", 1)
            if isinstance(res0, aa.PRO):
                self.apply_patch("NEUTRAL-NTERM", res0)
            elif neutraln:
                self.apply_patch("NEUTRAL-NTERM", res0)
            else:
                self.apply_patch("NTERM", res0)
        elif isinstance(res0, na.Nucleic):
            res0.set("is5term", 1)
            self.apply_patch("5TERM", res0)
        # Set the C-Terminus/ 3' Terminus
        reslast = chain.residues[-1]
        if isinstance(reslast, aa.Amino):
            reslast.set("is_c_terminus", 1)
            if neutralc:
                self.apply_patch("NEUTRAL-CTERM", reslast)
            else:
                self.apply_patch("CTERM", reslast)
        elif isinstance(reslast, na.Nucleic):
            reslast.set("is3term", 1)
            self.apply_patch("3TERM", reslast)
        else:
            for i in range(len(chain.residues)):
                resthis = chain.residues[-1 - i]
                if isinstance(resthis, aa.Amino):
                    resthis.set("is_c_terminus", 1)
                    if neutralc:
                        self.apply_patch("NEUTRAL-CTERM", resthis)
                    else:
                        self.apply_patch("CTERM", resthis)
                    break
                elif resthis.name in ["NH2", "NME"]:
                    break
                elif isinstance(resthis, na.Nucleic):
                    resthis.set("is3term", 1)
                    self.apply_patch("3TERM", resthis)
                    break
    def set_termini(self, neutraln=False, neutralc=False):
        """ Set the termini for the protein. First set all known termini by looking at the ends of
        the chain. Then examine each residue, looking for internal chain breaks. """
        self.write("Setting the termini... \n")
        # First assign the known termini
        for chain in self.protein.getChains():
            self.assign_termini(chain, neutraln, neutralc)
        # Now determine if there are any hidden chains
        letters = string.ascii_letters
        c = 0
        while c < len(self.protein.getChains()):
            chain = self.protein.chains[c]
            reslist = []
            origlist = []
            # origlist holds the original residue list for the chain
            for residue in chain.get_residues():
                origlist.append(residue)
            for residue in origlist:
                reslist.append(residue)
                # Look for ending termini
                fixflag = 0
                if isinstance(residue, aa.Amino):
                    if residue.has_atom("OXT") and not residue.is_c_terminus:
                        fixflag = 1
                elif isinstance(residue, na.Nucleic):
                    if (residue.has_atom("H3T") or residue.name.endswith("3")) \
                    and not residue.is3term:
                        fixflag = 1
                if fixflag:
                    # Get an available chain ID
                    chainid = letters[0]
                    id_ = 0
                    id_length = 1
                    while chainid in self.protein.chainmap:
                        id_ += 1
                        if id_ >= len(letters):
                            id_length += 1
                            id_ = 0
                        chainid = letters[id_] * id_length
                    if id_length > 1:
                        message = 'Warning: Reusing chain id: ' + chainid[0] + '\n'
                        self.write(message)
                    # Make a new chain with these residues
                    newchain = Chain(chainid[0])
                    self.protein.chainmap[chainid] = newchain
                    self.protein.chains.insert(c, newchain)
                    for res in reslist:
                        newchain.add_residue(res)
                        chain.residues.remove(res)
                        res.set_chain_id(chainid[0])
                    self.assign_termini(chain, neutraln, neutralc)
                    self.assign_termini(newchain, neutraln, neutralc)
                    reslist = []
                    c += 1
            c += 1
        # Update the final chain's chain_id if it is "" unless it's all water
        if "" in self.protein.chainmap:
            notwat = 0
            for res in chain.residues:
                if not isinstance(res, mols.WAT):
                    notwat = 1
                    break
            if notwat == 0:
                self.write("Done.\n")
                return
            chain = self.protein.chainmap[""]
            chainid = letters[0]
            id_ = 0
            id_length = 1
            while chainid in self.protein.chainmap:
                id_ += 1
                if id_ >= len(letters):
                    id_length += 1
                    id_ = 0
                chainid = letters[id_] * id_length
            if id_length > 1:
                message = 'Warning: Reusing chain id: ' + chainid[0] + '\n'
                self.write(message)
            # Use the new chain_id
            self.protein.chainmap[chainid] = chain
            del self.protein.chainmap[""]
            for res in chain.residues:
                res.set_chain_id(chainid[0])
        self.write("Done.\n")
    def find_missing_heavy(self):
        """ Repair residues that contain missing heavy (non-Hydrogen) atoms """
        self.write("Checking for missing heavy atoms... \n")
        misscount = 0
        heavycount = 0
        for residue in self.protein.get_residues():
            if not isinstance(residue, (aa.Amino, na.Nucleic)):
                continue
            # Check for Missing Heavy Atoms
            for refatomname in residue.reference.map:
                if refatomname.startswith("H"):
                    continue
                if refatomname in ["N+1", "C-1"]:
                    continue
                if refatomname in ["O1P", "O2P"]:
                    if residue.has_atom("OP1") and residue.has_atom("OP2"):
                        continue
                heavycount += 1
                if not residue.has_atom(refatomname):
                    self.write("Missing %s in %s\n" % \
                               (refatomname, residue), 1)
                    misscount += 1
                    residue.addMissing(refatomname)
            # Check for Extra Atoms
            atomlist = []
            for atom in residue.get("atoms"):
                atomlist.append(atom)
            for atom in atomlist:
                atomname = atom.get("name")
                if atomname in ["OP1", "OP2"] and residue.reference.has_atom("O1P") \
                    and residue.reference.has_atom("O2P"): continue
                if not residue.reference.has_atom(atomname):
                    self.write("Extra atom %s in %s! - " % \
                               (atomname, residue), 1)
                    residue.remove_atom(atomname)
                    self.write("Deleted this atom.\n")
        if heavycount == 0:
            raise PDBInputError("No heavy atoms found. You may also see this message if PDB2PQR \
does not have parameters for any residue in your protein.")
        misspct = 100.0 * float(misscount) / heavycount
        if misspct > REPAIR_LIMIT:
            error = "This PDB file is missing too many (%i out of " % misscount
            error += "%i, %.2f%%) heavy atoms to accurately repair the file.  " % \
                     (heavycount, misspct)
            error += "The current repair limit is set at %i%%. " % REPAIR_LIMIT
            error += "You may also see this message if PDB2PQR does not have parameters for enough \
residues in your protein."
            raise PDBInputError(error)
        elif misscount > 0:
            self.write("Missing %i out of %i heavy atoms (%.2f percent) - " % \
                       (misscount, heavycount, misspct))
            self.write("Will attempt to repair.\n")
            self.repair_heavy()
        else:
            self.write("No heavy atoms found missing - Done.\n")

    @staticmethod
    def rebuild_tetrahedral(residue, atomname):
        """ Rebuild a tetrahedral hydrogen group.  This is necessary due to the shortcomings of the
        quatfit routine - given a tetrahedral geometry and two existing hydrogens, the quatfit
        routines have two potential solutions.  This function uses basic tetrahedral geometry to fix
        this issue.
        Parameters
            residue:  The residue in question (residue)
            atomname: The atomname to add (string)
        Returns
            1 if successful, 0 otherwise """
        hcount = 0
        nextatomname = None
        atomref = residue.reference.map.get(atomname)
        if atomref is None:
            return False
        bondname = atomref.bonds[0]
        # Return if the bonded atom does not exist
        if not residue.has_atom(bondname):
            return False
        # This group is tetrahedral if bondatom has 4 bonds, 3 of which are hydrogens
        for bond in residue.reference.map[bondname].bonds:
            if bond.startswith("H"):
                hcount += 1
            elif bond != 'C-1' and bond != 'N+1':
                nextatomname = bond
        # Check if this is a tetrahedral group
        if hcount != 3 or nextatomname == None:
            return False
        # Now rebuild according to the tetrahedral geometry
        bondatom = residue.get_atom(bondname)
        nextatom = residue.get_atom(nextatomname)
        numbonds = len(bondatom.bonds)
        if numbonds == 1:
            # Place according to two atoms
            coords = [bondatom.getCoords(), nextatom.getCoords()]
            refcoords = [residue.reference.map[bondname].getCoords(), \
                         residue.reference.map[nextatomname].getCoords()]
            refatomcoords = atomref.getCoords()
            newcoords = find_coordinates(2, coords, refcoords, refatomcoords)
            residue.create_atom(atomname, newcoords)
            # For LEU and ILE residues only: make sure the Hydrogens are in staggered conformation
            # instead of eclipsed.
            if isinstance(residue, aa.LEU):
                hcoords = newcoords
                cbatom = residue.get_atom('CB')
                ang = get_dihedral(cbatom.getCoords(), nextatom.getCoords(), bondatom.getCoords(),
                                   hcoords)
                diffangle = 60 - ang
                residue.rotateTetrahedral(nextatom, bondatom, diffangle)
            elif isinstance(residue, aa.ILE):
                hcoords = newcoords
                cg1atom = residue.get_atom('CG1')
                cbatom = residue.get_atom('CB')
                if bondatom.name == 'CD1':
                    ang = get_dihedral(cbatom.getCoords(), nextatom.getCoords(),
                                       bondatom.getCoords(), hcoords)
                elif bondatom.name == 'CG2':
                    ang = get_dihedral(cg1atom.getCoords(), nextatom.getCoords(),
                                       bondatom.getCoords(), hcoords)
                else:
                    ang = get_dihedral(cbatom.getCoords(), nextatom.getCoords(),
                                       bondatom.getCoords(), hcoords)
                diffangle = 60 - ang
                residue.rotateTetrahedral(nextatom, bondatom, diffangle)
            return 1
        elif numbonds == 2:
            # Get the single hydrogen coordinates
            hatom = None
            for bond in bondatom.reference.bonds:
                if residue.has_atom(bond) and bond.startswith("H"):
                    hatom = residue.get_atom(bond)
                    break
            # Use the existing hydrogen and rotate about the bond
            residue.rotateTetrahedral(nextatom, bondatom, 120)
            newcoords = hatom.getCoords()
            residue.rotateTetrahedral(nextatom, bondatom, -120)
            residue.create_atom(atomname, newcoords)
            return 1
        elif numbonds == 3:
            # Find the one spot the atom can be
            hatoms = []
            for bond in bondatom.reference.bonds:
                if residue.has_atom(bond) and bond.startswith("H"):
                    hatoms.append(residue.get_atom(bond))
            # If this is more than two something is wrong
            if len(hatoms) != 2:
                return 0
            # Use the existing hydrogen and rotate about the bond
            residue.rotateTetrahedral(nextatom, bondatom, 120)
            newcoords1 = hatoms[0].getCoords()
            residue.rotateTetrahedral(nextatom, bondatom, 120)
            newcoords2 = hatoms[0].getCoords()
            residue.rotateTetrahedral(nextatom, bondatom, 120)
            # Determine which one hatoms[1] is not in
            if distance(hatoms[1].getCoords(), newcoords1) > 0.1:
                residue.create_atom(atomname, newcoords1)
            else:
                residue.create_atom(atomname, newcoords2)
            return 1

    def add_hydrogens(self):
        """ Add the hydrogens to the protein.  This requires either the rebuild_tetrahedral function
        for tetrahedral geometries or the standard quatfit methods.  These methods use three nearby
        bonds to rebuild the atom; the closer the bonds, the more accurate the results.  As such the
        peptide bonds are used when available. """
        count = 0
        self.write("Adding hydrogens to the protein...\n")
        for residue in self.protein.get_residues():
            if not isinstance(residue, (aa.Amino, na.Nucleic)):
                continue
            for atomname in residue.reference.map:
                if not atomname.startswith("H"):
                    continue
                if residue.has_atom(atomname):
                    continue
                if isinstance(residue, aa.CYS) and residue.disulfide_bonded and atomname == "HG":
                    continue
                # If this hydrogen is part of a tetrahedral group, follow a different codepath
                if Routines.rebuild_tetrahedral(residue, atomname):
                    count += 1
                    continue
                # Otherwise use the standard quatfit methods
                coords = []
                refcoords = []
                refatomcoords = residue.reference.map[atomname].getCoords()
                bondlist = residue.reference.get_nearest_bonds(atomname)
                for bond in bondlist:
                    if bond == "N+1":
                        atom = residue.peptide_n
                    elif bond == "C-1":
                        atom = residue.peptide_c
                    else:
                        atom = residue.get_atom(bond)
                    if atom == None:
                        continue
                    # Get coordinates, reference coordinates
                    coords.append(atom.getCoords())
                    refcoords.append(residue.reference.map[bond].getCoords())
                    # Exit if we have enough atoms
                    if len(coords) == 3:
                        break
                if len(coords) == 3:
                    newcoords = find_coordinates(3, coords, refcoords, refatomcoords)
                    residue.create_atom(atomname, newcoords)
                    count += 1
                else:
                    self.write("Couldn't rebuild %s in %s!\n" % (atomname, residue), 1)
        self.write(" Added %i hydrogen atoms.\n" % count)
    def remove_hydrogen(self):
        """ Remove hydrogens """
        self.write("Stripping hydrogens from the protein...\n")
        for residue in self.protein.get_residues():
            if not isinstance(residue, (aa.Amino, na.Nucleic)):
                continue
            for atom in residue.atoms[:]:
                if atom.isHydrogen():
                    residue.remove_atom(atom.name)
    def repair_heavy(self):
        """ Repair all heavy atoms.  Unfortunately the first time we get to an atom we might not be
        able to rebuild it - it might depend on other atoms to be rebuild first (think side chains).
        As such a 'seenmap' is used to keep track of what we've already seen and subsequent attempts
        to rebuild the atom. """
        self.write("Rebuilding missing heavy atoms... \n")
        for residue in self.protein.get_residues():
            if not isinstance(residue, (aa.Amino, na.Nucleic)):
                continue
            missing = residue.get("missing")
            if missing == []:
                continue
            # Initialize some variables
            seenmap = {}
            nummissing = len(missing)
            while len(missing) > 0:
                coords = []
                refcoords = []
                atomname = missing.pop(0)
                refatomcoords = residue.reference.map[atomname].getCoords()
                bondlist = residue.reference.get_nearest_bonds(atomname)
                for bond in bondlist:
                    if bond == "N+1":
                        atom = residue.peptide_n
                    elif bond == "C-1":
                        atom = residue.peptide_c
                    else: atom = residue.get_atom(bond)
                    if atom == None:
                        continue
                    # Get coordinates, reference coordinates
                    coords.append(atom.getCoords())
                    refcoords.append(residue.reference.map[bond].getCoords())
                    # Exit if we have enough atoms
                    if len(coords) == 3:
                        break
                # We might need other atoms to be rebuilt first
                if len(coords) < 3:
                    try:
                        seenmap[atomname] += 1
                    except KeyError:
                        seenmap[atomname] = 1
                    missing.append(atomname)
                    if seenmap[atomname] > nummissing:
                        text = "Too few atoms present to reconstruct or cap residue %s in \
structure!\n" % (residue)
                        text += "This error is generally caused by missing backbone atoms in this \
protein;\nyou must use an external program to complete gaps in the protein backbone.\n"
                        text += "Heavy atoms missing from %s: " % (residue)
                        text += ' '.join(missing)
                        raise PDBInputError(text)
                else: # Rebuild the atom
                    newcoords = find_coordinates(3, coords, refcoords, refatomcoords)
                    residue.create_atom(atomname, newcoords)
                    self.write("Added %s to %s at coordinates" % (atomname, residue), 1)
                    self.write(" %.3f %.3f %.3f\n" % \
                           (newcoords[0], newcoords[1], newcoords[2]))
        self.write("Done.\n")
    def set_reference_distance(self):
        """ Set the distance to the CA atom in the residue. This is necessary for determining which
        atoms are allowed to move during rotations.  Uses the shortest_path algorithm found in
        utilities.py. """
        for residue in self.protein.get_residues():
            if not isinstance(residue, aa.Amino):
                continue
            # Initialize some variables
            atom_map = {}
            caatom = residue.get_atom("CA")
            if caatom == None:
                text = "Cannot set references to %s without CA atom!\n"
                raise PDBInputError(text)
            # Set up the linked map
            for atom in residue.get_atoms():
                atom_map[atom] = atom.bonds
            # Run the algorithm
            for atom in residue.get_atoms():
                if atom.isBackbone():
                    atom.refdistance = -1
                elif residue.is_c_terminus and atom.name == "HO":
                    # special case for HO in Cterm
                    atom.refdistance = 3
                elif residue.is_n_terminus and (atom.name == "H3" or atom.name == "H2"):
                    # special case for H2 or H3 in Nterm
                    atom.refdistance = 2
                else:
                    atom.refdistance = len(shortest_path(atom_map, atom, caatom)) - 1
    def get_bump_score(self, residue):
        """ Get an bump score (for steric clashes) for the current structure. """
        # Do some setup
        self.cells = Cells(CELL_SIZE)
        self.cells.assignCells(self.protein)
        self.calculate_dihedral_angles()
        self.set_donors_and_acceptors()
        self.update_internal_linked_bonds()
        self.set_reference_distance()
        bumpscore = 0.0
        if not isinstance(residue, aa.Amino):
            return 0.0
        # Initialize variables
        for atom in residue.get_atoms():
            atomname = atom.name
            if atomname[0] != "H":
                continue
            bumpscore = bumpscore + self.get_bump_score_atom(atom)
        return bumpscore
    def get_bump_score_atom(self, atom):
        """ Find nearby atoms for conflict-checking.  Uses neighboring cells to compare atoms rather
        than an all versus all O(n^2) algorithm, which saves a great deal of time.  There are
        several instances where we ignore potential conflicts; these include donor/acceptor pairs,
        atoms in the same residue, and bonded CYS bridges.

        Parameters
            atom:  Find nearby atoms to this atom (Atom)
        Returns
            bumpscore: a bump score sum((dist-cutoff)**20 for all near atoms

        Jens rewrote this function from find_nearby_atoms to be usable for detecting bumps for
        optimizable hydrogens """
        # Initialize some variables
        residue = atom.residue
        atom_size = BUMP_HYDROGEN_SIZE if atom.isHydrogen() else BUMP_HEAVY_SIZE
        # Get atoms from nearby cells
        closeatoms = self.cells.getNearCells(atom)
        # Loop through and see if any are within the cutoff
        bumpscore = 0.0
        for closeatom in closeatoms:
            closeresidue = closeatom.residue
            if closeresidue == residue and (closeatom in atom.bonds or atom in closeatom.bonds):
                continue
            if not isinstance(closeresidue, aa.Amino):
                continue
            if isinstance(residue, aa.CYS):
                if residue.disulfide_bonded_partner == closeatom:
                    continue
            # Also ignore if this is a donor/acceptor pair
            pair_ignored = False
            if atom.isHydrogen() and len(atom.bonds) != 0 and atom.bonds[0].hdonor \
            and closeatom.hacceptor:
                continue
            if closeatom.isHydrogen() and len(closeatom.bonds) != 0 and closeatom.bonds[0].hdonor \
            and atom.hacceptor:
                continue
            dist = distance(atom.getCoords(), closeatom.getCoords())
            other_size = BUMP_HYDROGEN_SIZE if closeatom.isHydrogen() else BUMP_HEAVY_SIZE
            cutoff = atom_size + other_size
            if dist < cutoff:
                bumpscore = bumpscore + 1000.0
                if pair_ignored:
                    self.write('This bump is a donor/acceptor pair.\n')
        self.write('BUMPSCORE ' + str(bumpscore) + '\n')
        return bumpscore
    def debump_protein(self):
        """ Make sure that none of the added atoms were rebuilt on top of existing atoms.  See each
        called function for more information. """
        self.write("Checking if we must debump any residues... \n")
        # Do some setup
        self.cells = Cells(CELL_SIZE)
        self.cells.assignCells(self.protein)
        self.calculate_dihedral_angles()
        self.set_donors_and_acceptors()
        self.update_internal_linked_bonds()
        self.set_reference_distance()
        # Determine which residues to debump
        for residue in self.protein.get_residues():
            if not isinstance(residue, aa.Amino):
                continue
            # Initialize variables
            conflictnames = self.find_residue_conflicts(residue, True)
            if not conflictnames:
                continue
            # Otherwise debump the residue
            self.write("Starting to debump %s...\n" % residue, 1)
            self.write("Debumping cutoffs: %2.1f for heavy-heavy, %2.1f for hydrogen-heavy, and \
%2.1f for hydrogen-hydrogen.\n" % (BUMP_HEAVY_SIZE*2, BUMP_HYDROGEN_SIZE+BUMP_HEAVY_SIZE,
                                   BUMP_HYDROGEN_SIZE*2), 1)
            if self.debump_residue(residue, conflictnames):
                self.write("Debumping Successful!\n\n", 1)
            else:
                text = "WARNING: Unable to debump %s\n" % residue
                self.write("********\n%s********\n\n" % text)
                self.warnings.append(text)
        self.write("Done.\n")

    def find_residue_conflicts(self, residue, write_conflict_info=False):
        """ Find residues with conflicts"""
        conflictnames = []
        for atom in residue.get_atoms():
            atomname = atom.name
            if not atom.added:
                continue
            if atomname == "H":
                continue
            if atom.optimizeable:
                continue
            nearatoms = self.find_nearby_atoms(atom)
            # If something is too close, we must debump the residue
            if nearatoms != {}:
                conflictnames.append(atomname)
                if write_conflict_info:
                    for repatom in nearatoms:
                        self.write("%s %s is too close to %s %s\n" % \
                                   (residue, atomname, repatom.residue, repatom.name), 1)
        return conflictnames
    def score_dihedral_angle(self, residue, anglenum):
        """ Score the 'goodness' of a dihedral angle. """
        score = 0
        atomnames = residue.reference.dihedrals[anglenum].split()
        pivot = atomnames[2]
        moveablenames = self.get_moveable_names(residue, pivot)
        for name in moveablenames:
            nearatoms = self.find_nearby_atoms(residue.get_atom(name))
            for v in nearatoms.values():
                score += v
        return score
    def debump_residue(self, residue, conflictnames):
        """ Debump a specific residue.  Only should be called if the residue has been detected to
        have a steric clash conflict. If called, try to rotate about dihedral angles to resolve the
        conflict.

        Parameters
            residue:  The residue in question
            conflictnames:  A list of atomnames that were rebuilt too close to other atoms
        Returns
            True if successful, False otherwise """
        # Initialize some variables
        anglenum = -1
        current_conflict_names = conflictnames
        # Try (up to 10 times) to find a workable solution
        for _ in range(ANGLE_TEST_COUNT):
            anglenum = self.pick_dihedral_angle(residue, current_conflict_names, anglenum)
            if anglenum == -1:
                return False
            self.write("Using dihedral angle number %i to debump the residue.\n" % anglenum, 1)
            bestscore = self.score_dihedral_angle(residue, anglenum)
            found_improvement = False
            bestangle = original_angle = residue.dihedrals[anglenum]
            # Skip the first angle as it's already known.
            for i in range(1, ANGLE_STEPS):
                newangle = original_angle + (ANGLE_STEP_SIZE * i)
                self.set_dihedral_angle(residue, anglenum, newangle)
                # Check for conflicts
                score = self.score_dihedral_angle(residue, anglenum)
                if score == 0:
                    if not self.find_residue_conflicts(residue):
                        self.write("No conflicts found at angle "+repr(newangle)+"\n", 1)
                        return True
                    else:
                        bestangle = newangle
                        found_improvement = True
                        break
                # Set the best angle
                elif score < bestscore:
                    diff = abs(bestscore - score)
                    # Don't update if it's effectively a tie
                    if diff > EPSILON:
                        bestscore = score
                        bestangle = newangle
                        found_improvement = True
            self.set_dihedral_angle(residue, anglenum, bestangle)
            current_conflict_names = self.find_residue_conflicts(residue)
            if found_improvement:
                self.write("Best score of " + repr(bestscore) + " at angle " + repr(bestangle) +
                           ". New conflict set: ", 1)
                self.write(str(current_conflict_names)+"\n", 1)
            else:
                self.write("No improvement found for this dihedral angle.\n", 1)
        # If we're here, debumping was unsuccessful
        return False
    def calculate_dihedral_angles(self):
        """ Calculate the dihedral angle for every residue within the protein """
        for residue in self.protein.get_residues():
            if not isinstance(residue, aa.Amino):
                continue
            residue.dihedrals = []
            refangles = residue.reference.dihedrals
            for di in refangles:
                coords = []
                atoms = di.split()
                for i in range(4):
                    atomname = atoms[i]
                    if residue.has_atom(atomname):
                        coords.append(residue.get_atom(atomname).getCoords())
                if len(coords) == 4:
                    angle = get_dihedral(coords[0], coords[1], coords[2], coords[3])
                else:
                    angle = None
                residue.add_dihedral_angle(angle)
    def get_closest_atom(self, atom):
        """ Get the closest atom that does not form a donor/acceptor pair. Used to detect potential
        conflicts.

        NOTE:  Cells must be set before using this function.

        Parameters
            atom:  The atom in question (Atom)
        Returns
            bestatom:  The closest atom to the input atom that does not satisfy a donor/acceptor
            pair. """
        # Initialize some variables
        # TODO - The use of 999.99 drives me nuts.
        # TODO - Couldn't these be replaced with numpy.na; something easier to interpret at the end?
        bestdist = 999.99
        bestwatdist = 999.99
        bestatom = None
        bestwatatom = None
        residue = atom.residue
        # Get atoms from nearby cells
        closeatoms = self.cells.getNearCells(atom)
        # Loop through and see which is the closest
        for closeatom in closeatoms:
            closeresidue = closeatom.residue
            if closeresidue == residue:
                continue
            if not isinstance(closeresidue, (aa.Amino, mols.WAT)):
                continue
            if isinstance(residue, aa.CYS):
                if residue.disulfide_bonded_partner == closeatom:
                    continue
            # Also ignore if this is a donor/acceptor pair
            if atom.isHydrogen() and atom.bonds[0].hdonor and closeatom.hacceptor:
                continue
            if closeatom.isHydrogen() and closeatom.bonds[0].hdonor and atom.hacceptor:
                continue
            dist = distance(atom.getCoords(), closeatom.getCoords())
            if isinstance(closeresidue, mols.WAT):
                if dist < bestwatdist:
                    bestwatdist = dist
                    bestwatatom = closeatom
            else:
                if dist < bestdist:
                    bestdist = dist
                    bestatom = closeatom
        if bestdist > bestwatdist:
            txt = "Warning: %s in %s skipped when optimizing %s in %s\n" % (bestwatatom.name,
                                                                            bestwatatom.residue,
                                                                            atom.name, residue)
            if txt not in self.warnings:
                self.warnings.append(txt)
        return bestatom
    def find_nearby_atoms(self, atom):
        """ Find nearby atoms for conflict-checking.  Uses neighboring cells to compare atoms rather
        than an all versus all O(n^2) algorithm, which saves a great deal of time.  There are
        several instances where we ignore potential conflicts; these include donor/acceptor pairs,
        atoms in the same residue, and bonded CYS bridges.

        Parameters
            atom:  Find nearby atoms to this atom (Atom)
        Returns
            nearatoms:  A dictionary of <Atom too close> to <amount of overlap for that atom>.
        """
        # Initialize some variables
        nearatoms = {}
        residue = atom.residue
        atom_size = BUMP_HYDROGEN_SIZE if atom.isHydrogen() else BUMP_HEAVY_SIZE
        # Get atoms from nearby cells
        closeatoms = self.cells.getNearCells(atom)
        # Loop through and see if any are within the cutoff
        for closeatom in closeatoms:
            closeresidue = closeatom.residue
            if closeresidue == residue and (closeatom in atom.bonds or atom in closeatom.bonds):
                continue
            if not isinstance(closeresidue, (aa.Amino, mols.WAT)):
                continue
            if isinstance(residue, aa.CYS) and residue.disulfide_bonded_partner == closeatom:
                continue
            # Also ignore if this is a donor/acceptor pair
            if atom.isHydrogen() and (len(atom.bonds) != 0) \
            and atom.bonds[0].hdonor and closeatom.hacceptor:
                continue
            if closeatom.isHydrogen() and (len(closeatom.bonds) != 0) \
            and closeatom.bonds[0].hdonor and atom.hacceptor:
                continue
            dist = distance(atom.getCoords(), closeatom.getCoords())
            other_size = BUMP_HYDROGEN_SIZE if closeatom.isHydrogen() else BUMP_HEAVY_SIZE
            cutoff = atom_size + other_size
            if dist < cutoff:
                nearatoms[closeatom] = cutoff - dist
        return nearatoms

    def pick_dihedral_angle(self, residue, conflictnames, oldnum=None):
        """ Choose an angle number to use in debumping. Instead of simply picking a random chi
        angle, this function uses a more intelligent method to improve efficiency. The algorithm
        uses the names of the conflicting atoms within the residue to determine which angle number
        has the best chance of fixing the problem(s). The method also insures that the same chi
        angle will not be run twice in a row.

        Parameters
            residue: The residue that is being debumped (Residue)
            conflictnames: A list of atom names that are currently conflicts (list)
            oldnum: The old dihedral angle number (int)
        Returns the new dihedral angle number (int) """
        bestnum = -1
        best = 0
        i_list = range(len(residue.dihedrals))
        #Make sure our testing is done round robin.
        if oldnum is not None and oldnum >= 0 and len(i_list) > 0:
            del i_list[oldnum]
            test_dihedral_indices = i_list[oldnum:] + i_list[:oldnum]
        else:
            test_dihedral_indices = i_list
        for i in test_dihedral_indices:
            if i == oldnum:
                continue
            if residue.dihedrals[i] is None:
                continue
            score = 0
            atomnames = residue.reference.dihedrals[i].split()
            pivot = atomnames[2]
            moveablenames = self.get_moveable_names(residue, pivot)
            # If this pivot only moves the conflict atoms, pick it
            if conflictnames == moveablenames:
                return i
            # Otherwise find the pivot with the most matches
            for name in conflictnames:
                if name in moveablenames:
                    score += 1
                    if score > best:
                        best = score
                        bestnum = i
        # Return the best angle.  If none were found, return -1.
        return bestnum

    def set_dihedral_angle(self, residue, anglenum, angle):
        """ Rotate a residue about a given angle. Uses the quatfit methods to perform the matrix
        mathematics.
        Parameters
            residue: The residue to rotate
            anglenum: The number of the angle to rotate as listed in residue.dihedrals
            angle: The desired angle. """
        coordlist = []
        initcoords = []
        movecoords = []
        pivot = ""
        oldangle = residue.dihedrals[anglenum]
        diff = angle - oldangle
        atomnames = residue.reference.dihedrals[anglenum].split()
        pivot = atomnames[2]
        for atomname in atomnames:
            if residue.has_atom(atomname):
                coordlist.append(residue.get_atom(atomname).getCoords())
            else:
                raise PDBInputError("Error occurred while trying to debump!")
        initcoords = subtract(coordlist[2], coordlist[1])
        moveablenames = self.get_moveable_names(residue, pivot)
        for name in moveablenames:
            atom = residue.get_atom(name)
            movecoords.append(subtract(atom.getCoords(), coordlist[1]))
        newcoords = qchichange(initcoords, movecoords, diff)
        for i in range(len(moveablenames)):
            atom = residue.get_atom(moveablenames[i])
            self.cells.removeCell(atom)
            x = (newcoords[i][0] + coordlist[1][0])
            y = (newcoords[i][1] + coordlist[1][1])
            z = (newcoords[i][2] + coordlist[1][2])
            atom.set("x", x)
            atom.set("y", y)
            atom.set("z", z)
            self.cells.addCell(atom)
        # Set the new angle
        coordlist = []
        for atomname in atomnames:
            if residue.has_atom(atomname):
                coordlist.append(residue.get_atom(atomname).getCoords())
            else:
                raise PDBInputError("Error occurred while trying to debump!")
        di = get_dihedral(coordlist[0], coordlist[1], coordlist[2], coordlist[3])
        residue.dihedrals[anglenum] = di

    def get_moveable_names(self, residue, pivot):
        """ Return all atomnames that are further away than the pivot atom.
        Parameters
            residue:  The residue to use
            pivot:    The pivot atomname """
        movenames = []
        refdist = residue.get_atom(pivot).refdistance
        for atom in residue.get_atoms():
            if atom.refdistance > refdist:
                movenames.append(atom.name)
        return movenames

    def set_donors_and_acceptors(self):
        """ Set the donors and acceptors within the protein """
        for residue in self.protein.get_residues():
            residue.set_donors_and_acceptors()

    def run_pdb2pka(self, ph, ff, pdblist, ligand, verbose, pdb2pka_params):
        """ Run PDB2PKA """
        if ff.lower() != 'parse':
            PDB2PKAError('PDB2PKA can only be run with the PARSE force field.')
        self.write("Running PDB2PKA and applying at pH %.2f... \n" % ph)
        import pka
        from pdb2pka import pka_routines
        init_params = pdb2pka_params.copy()
        init_params.pop('pairene')
        init_params.pop('clean_output')
        results = pka.pre_init(original_pdb_list=pdblist, ff=ff, verbose=verbose, ligand=ligand,
                               **init_params)
        output_dir, protein, routines, forcefield, apbs_setup, \
        _, maps, sd = results
        my_pka_routines = pka_routines.pKaRoutines(protein, routines, forcefield, apbs_setup,
                                                   output_dir, maps, sd,
                                                   restart=pdb2pka_params.get('clean_output'),
                                                   pairene=pdb2pka_params.get('pairene'))
        print('Doing full pKa calculation')
        my_pka_routines.runpKa()
        pdb2pka_warnings = my_pka_routines.warnings[:]
        self.warnings.extend(pdb2pka_warnings)
        residue_ph = {}
        for pka_residue_tuple, calc_ph in my_pka_routines.ph_at_0_5.iteritems():
            tit_type, chain_id, number_str = pka_residue_tuple
            if tit_type == 'NTR':
                tit_type = 'N+'
            elif tit_type == 'CTR':
                tit_type = 'C-'
            key = ' '.join([tit_type, number_str, chain_id])
            residue_ph[key] = calc_ph
        pformat(residue_ph)
        self.apply_pka_values(ff, ph, residue_ph)
        self.write('Finished running PDB2PKA.\n')

    def run_propka(self, ph, ff, rootname, outname, options):
        """ Run PROPKA on the current protein, setting protonation states to the correct values
        Parameters
           ph:  The desired pH of the system
           ff:  The forcefield name to be used
           outname: The name of the PQR outfile """
        self.write("Running PROPKA and applying at pH %.2f... \n" % ph)
        from propka30.Source.protein import Protein as pkaProtein
        from propka30.Source.pdb import read_pdb as pkaReadPDB
        from propka30.Source.lib import residueList, setVerbose
        setVerbose(options.verbose)
        # Initialize some variables
        linelen = 70
        pkadic = {}
        # Reorder the atoms in each residue to start with N
        for residue in self.protein.get_residues():
            residue.reorder()
        # Make a string with all non-hydrogen atoms
        h_free_protein_file = StringIO()
        for atom in self.protein.get_atoms():
            if not atom.isHydrogen():
                atomtxt = atom.getPDBString()
                atomtxt = atomtxt[:linelen]
                h_free_protein_file.write(atomtxt)
                h_free_protein_file.write('\n')
        h_free_protein_file.seek(0)
        # Run PropKa
        atoms = pkaReadPDB('', file=h_free_protein_file)
        # creating protein object
        my_pka_protein = pkaProtein(atoms=atoms, name=rootname, options=options)
        # calculating pKa values for ionizable residues
        my_pka_protein.calculatePKA(options=options)
        # printing pka file
        my_pka_protein.writePKA(options=options, filename=outname)
        # Parse the results
        # This is the method used to generate the summary in the first place.
        residue_list = residueList("propka1")
        for chain in my_pka_protein.chains:
            for residue_type in residue_list:
                for residue in chain.residues:
                    if residue.res_name == residue_type:
                        #Strip out the extra space after C- or N+
                        key = '%s %s %s' % (residue.res_name.strip(), residue.resNumb,
                                            residue.chain_id)
                        key = key.strip()
                        pkadic[key] = residue.pKa_pro
        if len(pkadic) == 0:
            return
        # Now apply each pka to the appropriate residue
        self.apply_pka_values(ff, ph, pkadic)
        self.write("Done.\n")

    def apply_pka_values(self, ff, ph, pkadic):
        """ Apply pKa values to determine titration state at given pH """
        # TODO - This seems like a flaky way to do that since the titration curves (when available
        # have more information than the individual pKa values
        self.write('Applying pKa values at a pH of %.2f:\n' % ph)
        formatted_pkadict = pformat(pkadic)
        self.write(formatted_pkadict+'\n\n')
        warnings = []
        for residue in self.protein.get_residues():
            if not isinstance(residue, aa.Amino):
                continue
            resname = residue.name
            resnum = residue.res_seq
            chain_id = residue.chain_id
            if residue.is_n_terminus:
                key = "N+ %i %s" % (resnum, chain_id)
                key = key.strip()
                if key in pkadic:
                    value = pkadic[key]
                    del pkadic[key]
                    if ph >= value:
                        if ff in ["amber", "charmm", "tyl06", "peoepb", "swanson"]:
                            warn = ("N-terminal %s" % key, "neutral")
                            warnings.append(warn)
                        else:
                            self.apply_patch("NEUTRAL-NTERM", residue)
            if residue.is_c_terminus:
                key = "C- %i %s" % (resnum, chain_id)
                key = key.strip()
                if key in pkadic:
                    value = pkadic[key]
                    del pkadic[key]
                    if ph < value:
                        if ff in ["amber", "charmm", "tyl06", "peoepb", "swanson"]:
                            warn = ("C-terminal %s" % key, "neutral")
                            warnings.append(warn)
                        else:
                            self.apply_patch("NEUTRAL-CTERM", residue)
            key = "%s %i %s" % (resname, resnum, chain_id)
            key = key.strip()
            if key in pkadic:
                value = pkadic[key]
                del pkadic[key]
                if resname == "ARG" and ph >= value:
                    if ff == "parse":
                        self.apply_patch("AR0", residue)
                        txt = "WARNING: Neutral arginines are very rare. Please double\n"
                        self.warnings.append(txt)
                        self.write(txt)
                        txt = "         check your system and caculation setup.\n"
                        self.warnings.append(txt)
                        self.write(txt)
                    else:
                        warn = (key, "neutral")
                        warnings.append(warn)
                elif resname == "ASP" and ph < value:
                    if residue.is_c_terminus and ff in ["amber", "tyl06", "swanson"]:
                        warn = (key, "Protonated at C-Terminal")
                        warnings.append(warn)
                    elif residue.is_n_terminus and ff in ["amber", "tyl06", "swanson"]:
                        warn = (key, "Protonated at N-Terminal")
                        warnings.append(warn)
                    else:
                        self.apply_patch("ASH", residue)
                elif resname == "CYS" and ph >= value:
                    if ff == "charmm":
                        warn = (key, "negative")
                        warnings.append(warn)
                    else:
                        self.apply_patch("CYM", residue)
                elif resname == "GLU" and ph < value:
                    if residue.is_c_terminus and ff in ["amber", "tyl06", "swanson"]:
                        warn = (key, "Protonated at C-Terminal")
                        warnings.append(warn)
                    elif residue.is_n_terminus and ff in ["amber", "tyl06", "swanson"]:
                        warn = (key, "Protonated at N-Terminal")
                        warnings.append(warn)
                    else:
                        self.apply_patch("GLH", residue)
                elif resname == "HIS" and ph < value:
                    self.apply_patch("HIP", residue)
                elif resname == "LYS" and ph >= value:
                    if ff == "charmm":
                        warn = (key, "neutral")
                        warnings.append(warn)
                    elif ff in ["amber", "tyl06", "swanson"] and residue.get("is_c_terminus"):
                        warn = (key, "neutral at C-Terminal")
                        warnings.append(warn)
                    elif ff == "tyl06" and residue.get("is_n_terminus"):
                        warn = (key, "neutral at N-Terminal")
                        warnings.append(warn)
                    else:
                        self.apply_patch("LYN", residue)
                elif resname == "TYR" and ph >= value:
                    if ff in ["charmm", "amber", "tyl06", "peoepb", "swanson"]:
                        warn = (key, "negative")
                        warnings.append(warn)
                    else:
                        self.apply_patch("TYM", residue)
        if len(warnings) > 0:
            init = "WARNING: PDB2PKA determined the following residues to be\n"
            self.warnings.append(init)
            self.write(init)
            init = "         in a protonation state not supported by the\n"
            self.warnings.append(init)
            self.write(init)
            init = "         %s forcefield!\n" % ff
            self.warnings.append(init)
            self.write(init)
            init = "         All were reset to their standard pH 7.0 state.\n"
            self.warnings.append(init)
            self.write(init)
            self.warnings.append("\n")
            self.write('\n')
            for warn in warnings:
                text = "             %s (%s)\n" % (warn[0], warn[1])
                self.warnings.append(text)
                self.write(text)
            self.warnings.append("\n")
            self.write('\n')
        if len(pkadic) > 0:
            warn = "         PDB2PQR could not identify the following residues\n"
            self.warnings.append(warn)
            self.write(warn)
            warn = "         and residue numbers as returned by PROPKA or PDB2PKA:\n"
            self.warnings.append(warn)
            self.warnings.append("\n")
            self.write(warn)
            self.write('\n')
            for item in pkadic:
                text = "             %s\n" % item
                self.warnings.append(text)
                self.write(text)
            self.warnings.append("\n")
            self.write('\n')

class Cells:
    """ The cells object provides a better way to search for nearby atoms. A pure all versus all
    search is O(n^2) - for every atom, every other atom must be searched.  This is rather
    inefficient, especially for large proteins where cells may be tens of angstroms apart.  The cell
    class breaks down the xyz protein space into several 3-D cells of desired size - then by simply
    examining atoms that fall into the adjacent cells one can quickly find nearby cells. """
    # TODO - this is more general than everything else in the routines module; it should be moved
    def __init__(self, cellsize):
        """ Parameters
            cellsize:  The size of each cell (int) """
        self.cellmap = {}
        self.cellsize = cellsize

    def assignCells(self, protein):
        """ Place each atom in a virtual cell for easy neighbor comparison """
        for atom in protein.get_atoms():
            atom.cell = None
            self.addCell(atom)

    def addCell(self, atom):
        """ Add an atom to the cell
        Parameters
            atom:  The atom to add (atom) """
        size = self.cellsize
        x = atom.get("x")
        if x < 0:
            x = (int(x) - 1) / size * size
        else:
            x = int(x) / size * size
        y = atom.get("y")
        if y < 0:
            y = (int(y) - 1) / size * size
        else:
            y = int(y) / size * size
        z = atom.get("z")
        if z < 0:
            z = (int(z) - 1) / size * size
        else:
            z = int(z) / size * size
        key = (x, y, z)
        try:
            self.cellmap[key].append(atom)
        except KeyError:
            self.cellmap[key] = [atom]
        atom.set("cell", key)

    def removeCell(self, atom):
        """ Remove the atom from a cell
        Parameters
             atom:   The atom to add (atom) """
        oldcell = atom.get("cell")
        if oldcell == None:
            return
        atom.set("cell", None)
        self.cellmap[oldcell].remove(atom)

    def getNearCells(self, atom):
        """ Find all atoms in bordering cells to an atom
        Parameters
            atom:  The atom to use (atom)
            Returns a list of nearby atoms (list) """
        size = self.cellsize
        closeatoms = []
        cell = atom.get("cell")
        if cell == None:
            return closeatoms
        else:
            x = cell[0]
            y = cell[1]
            z = cell[2]
            for i in range(-1 * size, 2 * size, size):
                for j in range(-1 * size, 2 * size, size):
                    for k in range(-1 * size, 2 * size, size):
                        newkey = (x + i, y + j, z + k)
                        try:
                            newatoms = self.cellmap[newkey]
                            for atom2 in newatoms:
                                if atom == atom2:
                                    continue
                                closeatoms.append(atom2)
                        except KeyError:
                            pass
            return closeatoms
