""" Utilities for PDB2PQR Suite

    This module provides various utilities for the PDB2PQR suite to be imported into other Python
    scripts.

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

# TODO - Where do these numbers come from?
SMALL = 1.0e-7
DIHEDRAL = 57.2958

import math
import os
from os.path import splitext
import sys
import urllib.request, urllib.parse, urllib.error


class PROPKAoptions():
    """ Options for running PROPKA """
    def __init__(self, pH, verbose=False, reference='neutral'):
        """ Create a propka options object for running propka. """
        self.pH = pH
        self.reference = reference
        self.chains = None
        self.thermophiles = None
        self.alignment = None
        self.mutations = None
        self.verbose = verbose
        self.protonation = "old-school"
        self.window = (0.0, 14.0, 1.0)
        self.grid = (0.0, 14.0, 0.1)
        self.mutator = None
        self.mutator_options = None
        self.display_coupled_residues = None
        self.print_iterations = None
        self.version_label = "Nov30"
        # TODO - Make sure that PROPKA gets updated to 3.1 appropriately...
        from propka30.Source import lib
        lib.interpretMutator(self)
        lib.setDefaultAlignmentFiles(self)

def get_pqr_base_file_name(filename):
    """ Remove the stem """
    root, ext = splitext(filename)
    if ext.lower() == '.pqr':
        return root
    return filename

def sort_dict_by_value(inputdict):
    """ Sort a dictionary by its values """
    # TODO - go away or I will replace you with a lambda function!
    items = [(v, k) for k, v in list(inputdict.items())]
    items.sort()
    items.reverse()
    items = [k for v, k in items]
    return items

def shortest_path(graph, start, end, path=[]):
    """ Uses recursion to find the shortest path from one node to another in an unweighted graph.
    Adapted from http://www.python.org/doc/essays/graphs.html .
    Parameters:
    * graph: A mapping of the graph to analyze, of the form {0: [1,2], 1:[3,4], ...}. Each key has a
      list of edges.
    * start: The ID of the key to start the analysis from
    * end: The ID of the key to end the analysis
    * path: Optional argument used during the recursive step to keep the current path up to that
      point
    Returns a list of the shortest path (list) Returns None if start and end are not connected """
    # TODO - it would be better to replace this function with something from the networkx module
    path = path + [start]
    if start == end:
        return path
    if start not in graph:
        return None
    shortest = None
    for node in graph[start]:
        if node not in path:
            newpath = shortest_path(graph, node, end, path)
            if newpath:
                if not shortest or len(newpath) < len(shortest):
                    shortest = newpath
    return shortest

def analyze_connectivity(in_map, key):
    """ Analyze the connectivity of a given map using the key value.
    Parameters
        map:  The map to analyze (dict)
        key:  The key value (variable)
    Returns
        list: A list of connected values to the key (list) """
    key_list = []
    keys = [key]
    while len(keys) > 0:
        key = keys[0]
        if key not in key_list:
            key_list.append(key)
            # The following 4 lines are modified by Greg Cipriano as a bug fix
            if key in in_map:
                for value in in_map[key]:
                    if value not in key_list:
                        keys.append(value)
        keys.pop(keys.index(key))
    return list

def get_angle(coords1, coords2, coords3):
    """ Get the angle between three coordinates
    Parameters
        coords1:  The first coordinate set (atom)
        coords2:  The second (vertex) coordinate set (atom)
        coords3:  The third coordinate set (atom)
    Returns the angle between the atoms (float) """
    angle = 0.0
    v3m2 = subtract(coords3, coords2)
    v1m2 = subtract(coords1, coords2)
    norm1 = normalize(v3m2)
    norm2 = normalize(v1m2)
    dotted = dot(norm1, norm2)
    if dotted > 1.0: # If normalized, this is due to rounding error
        dotted = 1.0
    rad = abs(math.acos(dotted))
    angle = rad*180.0/math.pi
    if angle > 180.0:
        angle = 360.0 - angle
    return angle

def get_forcefield_file(name):
    """ Grab the forcefield file.  May or may not residue in the dat/ directory. """
    if name is None:
        return ''
    path = ""
    dirs = sys.path + ["dat"]
    if name in ["amber", "charmm", "parse", "tyl06", "peoepb", "swanson"]:
        name = name.upper()
    names = ["dat/%s.DAT" % name]
    names.append("%s.DAT" % name)
    names.append("%s.dat" % name)
    names.append("dat/%s" % name)
    names.append(name)
    for guess in names:
        if os.path.isfile(guess):
            return guess
        for path in dirs:
            testpath = "%s/%s" % (path, guess)
            if os.path.isfile(testpath):
                return testpath
    # If we get here return empty string
    return ""

def get_names_file(name):
    """ Grab the *.names file that contains the XML mapping.
    Parameters
        name:  The name of the forcefield (string)
    Returns the path to the file (string)"""
    if name is None:
        return ''
    path = ""
    dirs = sys.path + ["dat"]
    if name in ["amber", "charmm", "parse", "tyl06", "peoepb", "swanson"]:
        name = name.upper()
    names = ["dat/%s.names" % name]
    names.append("%s.names" % name)
    for guess in names:
        if os.path.isfile(guess):
            return guess
        for path in dirs:
            testpath = "%s/%s" % (path, guess)
            if os.path.isfile(testpath):
                return testpath
    # If we get here return empty string
    return ""

def get_data_file(name):
    """ Grab a data file. If the file cannot be found in the given directory, try the current system
    path.
    Parameters
        name:  The name of the file to get (string)
    Returns the path to the file (string) """
    path = ""
    if os.path.isfile(name):
        path = name
    for path in sys.path:
        testpath = "%s/%s" % (path, name)
        if os.path.isfile(testpath):
            path = testpath
    return path

def get_pdb_file(path):
    """ Obtain a PDB file.  First check the path given on the command line - if that file is not
    available, obtain the file from the PDB webserver at http://www.rcsb.org/pdb/ .
    Parameters
        path:  Name of PDB file to obtain (string)
    Returns file object containing PDB file (file object) """
    file = None
    if not os.path.isfile(path):
        urlpath = "http://www.rcsb.org/pdb/cgi/export.cgi/" + path + \
                  ".pdb?format=PDB&pdbId=" + path + "&compression=None"
        try:
            file = urllib.request.urlopen(urlpath)
            if file.getcode() != 200 or 'nosuchfile' in file.geturl():
                raise IOError
        except IOError:
            return None
    else:
        file = open(path, 'rU')
    return file

def distance(coords1, coords2):
    """ Calculate the distance between two coordinates
    Parameters
        coords1: Coordinates of form [x,y,z]
        coords2: Coordinates of form [x,y,z]
    Returns distance between the two coordinates (float) """
    # TODO - could be replaced with NumPy
    dist = 0.0
    dx = coords2[0] - coords1[0]
    dy = coords2[1] - coords1[1]
    dz = coords2[2] - coords1[2]
    dist = math.sqrt(dx*dx + dy*dy + dz*dz)
    return dist

def add(coords1, coords2):
    """ Add one 3-dimensional point to another
    Parameters
        coords1: coordinates of form [x,y,z]
        coords2: coordinates of form [x,y,z]
    Returns list of coordinates equal to coords2 + coords1 (list) """
    # TODO - could be replaced with NumPy
    x = coords1[0] + coords2[0]
    y = coords1[1] + coords2[1]
    z = coords1[2] + coords2[2]
    return [x, y, z]

def subtract(coords1, coords2):
    """ Subtract one 3-dimensional point from another
    Parameters
        coords1: coordinates of form [x,y,z]
        coords2: coordinates of form [x,y,z]
    Returns list of coordinates equal to coords1 - coords2 (list) """
    # TODO - could be replaced with NumPy
    x = coords1[0] - coords2[0]
    y = coords1[1] - coords2[1]
    z = coords1[2] - coords2[2]
    return [x, y, z]

def cross(coords1, coords2):
    """ Find the cross product of two 3-dimensional points
    Parameters
        coords1: coordinates of form [x,y,z]
        coords2: coordinates of form [x,y,z]
    Returns cross product coords2 and coords1 (list) """
    # TODO - could be replaced with NumPy
    x = coords1[1]*coords2[2] -  coords1[2]*coords2[1]
    y = coords1[2]*coords2[0] -  coords1[0]*coords2[2]
    z = coords1[0]*coords2[1] -  coords1[1]*coords2[0]
    return [x, y, z]

def dot(coords1, coords2):
    """ Find the dot product of two 3-dimensional points
    Parameters
        coords1: coordinates of form [x,y,z]
        coords2: coordinates of form [x,y,z]
    Returns dot product coords2 and coords1 (float) """
    # TODO - should be replaced with NumPy
    value = 0.0
    for i in range(3):
        value += coords1[i]*coords2[i]
    return value

def normalize(coords):
    """ Normalize a set of coordinates
    Parameters
        coords: coordinates of form [x,y,z]
    Returns normalized coordinates (list) """
    # TODO - could be replaced with NumPy
    dist = math.sqrt(pow(coords[0], 2) + pow(coords[1], 2) + pow(coords[2], 2))
    if dist > SMALL:
        nx = coords[0]/dist
        ny = coords[1]/dist
        nz = coords[2]/dist
        return [nx, ny, nz]
    else:
        return coords

def factorial(num):
    """ Returns the factorial of the given number n """
    # TODO - should be replaced with NumPy
    if num <= 1:
        return 1
    return num*factorial(num-1)

def get_dihedral(coords1, coords2, coords3, coords4):
    """ Calculate the angle using the four atoms
    Parameters
        coords1: First of four coordinates of form [x,y,z]
        coords2: Second of four
        coords3: Third of four
        coords4: Fourth of four
    Returns size of the angle (float) """
    value = 0.0
    v4m3 = subtract(coords4, coords3)
    v3m2 = subtract(coords3, coords2)
    v1m2 = subtract(coords1, coords2)
    v1232 = cross(v1m2, v3m2)
    v1232 = normalize(v1232)
    v4332 = cross(v4m3, v3m2)
    v4332 = normalize(v4332)
    scal = dot(v1232, v4332)
    if abs(scal + 1.0) < SMALL:
        value = 180.0
    elif abs(scal - 1.0) < SMALL:
        value = 0.0
    else:
        value = DIHEDRAL * math.acos(scal)
    chiral = dot(cross(v1232, v4332), v3m2)
    if chiral < 0:
        value = value * -1.0
    return value
