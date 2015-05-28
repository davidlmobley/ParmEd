#!/usr/bin/env python

import chemistry as chem
from simtk.openmm import *; from simtk.openmm.app import *; from simtk.unit import *

#Jason's original version - AMBER vs GROMACS with same (original) coordinates)
l = chem.amber.AmberParm('toluene_cyclohexane_10_500.prmtop', 'toluene_cyclohexane_10_500.inpcrd')
k = chem.load_file('toluene_cyclohexane_10_500_parmed.top')
[a.name for a in l.atoms] == [a.name for a in k.atoms]
sysl = l.createSystem()
sysk = k.createSystem()
conl = Context(sysl, VerletIntegrator(0.001))
conk = Context(sysk, VerletIntegrator(0.001))
conl.setPositions(l.positions)
conk.setPositions(l.positions)
print(chem.openmm.utils.energy_decomposition(l, conl))
print(chem.openmm.utils.energy_decomposition(k, conk))

#My new version - AMBER vs GROMACS with distinct (but same precision) coordinates generated from original coordinates
print("\nAMBER vs GROMACS with distinct (but same precision) coordinates:")
amberparm = chem.amber.AmberParm( 'toluene_cyclohexane_10_500_parmed.prmtop', 'toluene_cyclohexane_10_500_parmed.crd')
gromacsparm = chem.load_file( 'toluene_cyclohexane_10_500_parmed.top' )
gro = chem.gromacs.GromacsGroFile.parse( 'toluene_cyclohexane_10_500_parmed.gro' )
#Store coordinates, box
gromacsparm.box = gro.box
gromacsparm.positions = gro.positions
#Check atom names
[ a.name for a in amberparm.atoms ] == [ a.name for a in gromacsparm.atoms ]
#Create systems
sysa = amberparm.createSystem()
sysg = gromacsparm.createSystem()
cona = Context( sysa, VerletIntegrator(0.001))
cong = Context( sysg, VerletIntegrator(0.001))
cona.setPositions(amberparm.positions)
cong.setPositions(gromacsparm.positions)
print(chem.openmm.utils.energy_decomposition(amberparm, cona))
print(chem.openmm.utils.energy_decomposition(gromacsparm, cong))


#Could generalize by making a whole matrix of systems...
