"""
Contains classes for parsing GROMACS topology and parameter files
"""
import os as _os

#from chemistry.gromacs.gromacstop import GromacsTopologyFile
from chemistry.utils import which as _which

GROMACS_TOPDIR = None

if _os.getenv('GMXDATA') is not None and _os.path.isdir(
        _os.path.join(_os.getenv('GMXDATA'), 'top')):
    GROMACS_TOPDIR = _os.path.join(_os.getenv('GMXDATA'), 'top')
elif _os.getenv('GMXBIN') is not None and _os.path.isdir(
        _os.path.join(_os.getenv('GMXBIN'), '..', 'share', 'gromacs', 'top')):
    GROMACS_TOPDIR = _os.path.join(_os.getenv('GMXBIN'), '..', 'share',
                                  'gromacs', 'top')
else:
    for _testdir in ['/usr', '/usr/local', '/opt/local', '/opt']:
        if _os.path.isdir(_os.path.join(_testdir, 'share', 'gromacs')):
            GROMACS_TOPDIR = _os.path.join(_testdir, 'share', 'gromacs', 'top')
            break

if GROMACS_TOPDIR is None:
    if _which('gmx') is not None:
        GROMACS_TOPDIR = _os.path.join(_os.path.split(_which('gmx'))[0],
                                       '..', 'share', 'gromacs', 'top')

if GROMACS_TOPDIR is not None:
    # Regularize the include path
    GROMACS_TOPDIR = _os.path.realpath(GROMACS_TOPDIR)
else:
    # Use the default Gromacs installation path
    GROMACS_TOPDIR = '/usr/local/gromacs/share/gromacs/top'

try:
    del _testdir, _os, _which
except NameError: #_testdir is only defined under certain scenarios so it may not be defined here
    del _os, _which

from chemistry.gromacs.gromacstop import GromacsTopologyFile
from chemistry.gromacs.gromacsgro import GromacsGroFile
