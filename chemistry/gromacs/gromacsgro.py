"""
This module contains functionality relevant to loading and parsing GROMACS GRO
(coordinate) files and building a stripped-down Structure from it
"""
from __future__ import print_function, division, absolute_import

from chemistry.constants import TINY
from chemistry.exceptions import ParsingError
from chemistry.formats.registry import FileFormatType
from chemistry.geometry import (box_vectors_to_lengths_and_angles,
                                box_lengths_and_angles_to_vectors,
                                reduce_box_vectors)
from chemistry.structure import Structure
from chemistry.topologyobjects import Atom
from chemistry import unit as u
from chemistry.utils.io import genopen
from chemistry.utils.six import add_metaclass, string_types
from contextlib import closing
try:
    import numpy as np
except ImportError:
    np = None

@add_metaclass(FileFormatType)
class GromacsGroFile(object):
    """ Parses and writes Gromacs GRO files """
    #===================================================

    @staticmethod
    def id_format(filename):
        """ Identifies the file as a GROMACS GRO file

        Parameters
        ----------
        filename : str
            Name of the file to check if it is a Gromacs GRO file

        Returns
        -------
        is_fmt : bool
            If it is identified as a Gromacs GRO file, return True. False
            otherwise
        """
        with closing(genopen(filename)) as f:
            f.readline() # Title line
            try:
                int(f.readline().strip()) # number of atoms
            except ValueError:
                return False
            line = f.readline()
            try:
                int(line[:5])
                if not line[5:10].strip(): return False
                if not line[10:15].strip(): return False
                int(line[15:20])
                pdeci = [i for i, x in enumerate(line) if x == '.']
                ndeci = pdeci[1] - pdeci[0] - 5
                for i in range(1, 4):
                    wbeg = (pdeci[0]-4)+(5+ndeci)*(i-1)
                    wend = (pdeci[0]-4)+(5+ndeci)*i
                    float(line[wbeg:wend])
                i = 4
                wbeg = (pdeci[0]-4)+(5+ndeci)*(i-1)
                wend = (pdeci[0]-4)+(5+ndeci)*i
                if line[wbeg:wend].strip():
                    for i in range(4, 7):
                        wbeg = (pdeci[0]-4)+(5+ndeci)*(i-1)
                        wend = (pdeci[0]-4)+(5+ndeci)*i
                        float(line[wbeg:wend])
            except ValueError:
                return False
            return True

    #===================================================

    @staticmethod
    def parse(filename):
        """ Parses a Gromacs GRO file

        Parameters
        ----------
        filename : str or file-like
            Name of the file or the GRO file object

        Returns
        -------
        struct : :class:`Structure`
            The Structure instance instantiated with *just* residues and atoms
            populated (with coordinates)
        """
        struct = Structure()
        if isinstance(filename, string_types):
            fileobj = genopen(filename, 'r')
            own_handle = True
        else:
            fileobj = filename
            own_handle = False
        try:
            # Ignore the title line
            fileobj.readline()
            try:
                natom = int(fileobj.readline().strip())
            except ValueError:
                raise ParsingError('Could not parse %s as GRO file' % filename)
            for i, line in enumerate(fileobj):
                if i == natom: break
                try:
                    resnum = int(line[:5])
                    resname = line[5:10].strip()
                    atomname = line[10:15].strip()
                    atnum = int(line[15:20])
                    atom = Atom(name=atomname, number=atnum)
                    pdeci = [i for i, x in enumerate(line) if x == '.']
                    ndeci = pdeci[1] - pdeci[0] - 5
                    for i in range(1, 4):
                        wbeg = (pdeci[0]-4)+(5+ndeci)*(i-1)
                        wend = (pdeci[0]-4)+(5+ndeci)*i
                        if i == 1: atom.xx = float(line[wbeg:wend]) * 10
                        elif i == 2: atom.xy = float(line[wbeg:wend]) * 10
                        elif i == 3: atom.xz = float(line[wbeg:wend]) * 10
                    i = 4
                    wbeg = (pdeci[0]-4)+(5+ndeci)*(i-1)
                    wend = (pdeci[0]-4)+(5+ndeci)*i
                    if line[wbeg:wend].strip():
                        for i in range(4, 7):
                            wbeg = (pdeci[0]-4)+(5+ndeci)*(i-1)
                            wend = (pdeci[0]-4)+(5+ndeci)*i
                            if i == 4: atom.vx = float(line[wbeg:wend]) * 10
                            elif i == 5: atom.vy = float(line[wbeg:wend]) * 10
                            elif i == 5: atom.vz = float(line[wbeg:wend]) * 10
                except (ValueError, IndexError):
                    raise ParsingError('Could not parse the atom record of '
                                       'GRO file %s' % filename)
                struct.add_atom(atom, resname, resnum)
            # Get the box from the last line if it's present
            if line.strip():
                try:
                    box = [float(x) for x in line.split()]
                except ValueError:
                    raise ParsingError('Could not understand box line of GRO '
                                       'file %s' % filename)
                if len(box) == 3:
                    struct.box = [box[0]*10, box[1]*10, box[2]*10,
                                  90.0, 90.0, 90.0]
                elif len(box) == 9:
                    # Assume we have vectors
                    leng, ang = box_vectors_to_lengths_and_angles(
                                [box[0], box[3], box[4]]*u.nanometers,
                                [box[5], box[1], box[6]]*u.nanometers,
                                [box[7], box[8], box[2]]*u.nanometers)
                    a, b, c = leng.value_in_unit(u.angstroms)
                    alpha, beta, gamma = ang.value_in_unit(u.degrees)
                    struct.box = [a, b, c, alpha, beta, gamma]
                if np is not None:
                    struct.box = np.array(struct.box)
        finally:
            if own_handle:
                fileobj.close()

        return struct

    #===================================================

    @staticmethod
    def write(struct, dest, nobox=False):
        """ Write a Gromacs Topology File from a Structure

        Parameters
        ----------
        struct : :class:`Structure`
            The structure to write to a Gromacs GRO file (must have coordinates)
        dest : str or file-like
            The name of a file or a file object to write the Gromacs topology to
        nobox : bool, optional
            If the system does not have a periodic box defined, and this option
            is True, no box will be written. If False, the periodic box will be
            defined to enclose the solute with 0.5 nm clearance on all sides. If
            periodic box dimensions *are* defined, this variable has no effect.
        """
        if isinstance(dest, string_types):
            dest = genopen(dest, 'w')
            own_handle = True
        elif not hasattr(dest, 'write'):
            raise TypeError('dest must be a file name or file-like object')

        dest.write('GROningen MAchine for Chemical Simulation\n')
        dest.write('%5d\n' % len(struct.atoms))
        has_vels = all(hasattr(a, 'vx') for a in struct.atoms)
        for atom in struct.atoms:
            if has_vels:
                dest.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n' %
                           (atom.residue.idx+1, atom.residue.name, atom.name,
                               atom.idx+1, atom.xx/10, atom.xy/10, atom.xz/10,
                               atom.vx/10, atom.vy/10, atom.vz/10))
            else:
                dest.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %
                           (atom.residue.idx+1, atom.residue.name, atom.name,
                               atom.idx+1, atom.xx/10, atom.xy/10, atom.xz/10))
        # Box, in the weird format...
        if struct.box is not None:
            a, b, c = reduce_box_vectors(*box_lengths_and_angles_to_vectors(
                            *struct.box))
            if all([abs(x-90) < TINY for x in struct.box[3:]]):
                dest.write('%10.5f'*3 % (a[0]/10, b[1]/10, c[2]/10))
            else:
                dest.write('%10.5f'*9 % (a[0]/10, b[1]/10, c[2]/10, a[1]/10,
                           a[2]/10, b[0]/10, b[2]/10, c[0]/10, c[1]/10))
            dest.write('\n')
        elif not nobox and struct.atoms:
            # Find the extent of the molecule in all dimensions
            xdim = [struct.atoms[0].xx, struct.atoms[1].xx]
            ydim = [struct.atoms[0].xy, struct.atoms[1].xy]
            zdim = [struct.atoms[0].xz, struct.atoms[1].xz]
            for atom in struct.atoms:
                xdim[0] = min(xdim[0], atom.xx)
                xdim[1] = max(xdim[1], atom.xx)
                ydim[0] = min(ydim[0], atom.xy)
                ydim[1] = max(ydim[1], atom.xy)
                zdim[0] = min(zdim[0], atom.xz)
                zdim[1] = max(zdim[1], atom.xz)
            dest.write('%10.5f'*3 % (xdim[1]-xdim[0]+5, ydim[1]-ydim[0]+5,
                                     zdim[1]-zdim[0]+5))
            dest.write('\n')
        if own_handle:
            dest.close()
