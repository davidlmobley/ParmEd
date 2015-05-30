#!/usr/bin/env python

import chemistry as chem
from simtk.openmm import *; from simtk.openmm.app import *; from simtk.unit import *
import simtk.unit as u

#Jason's original version - AMBER vs GROMACS with same (original) coordinates)
l = chem.amber.AmberParm('toluene_cyclohexane_10_500.prmtop', 'toluene_cyclohexane_10_500.inpcrd')
k = chem.load_file('toluene_cyclohexane_10_500_parmed.top')
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
#Create systems
sysa = amberparm.createSystem()
sysg = gromacsparm.createSystem()
cona = Context( sysa, VerletIntegrator(0.001))
cong = Context( sysg, VerletIntegrator(0.001))
cona.setPositions(amberparm.positions)
cong.setPositions(gromacsparm.positions)
print(chem.openmm.utils.energy_decomposition(amberparm, cona))
print(chem.openmm.utils.energy_decomposition(gromacsparm, cong))

#Same as first case but with PME/PBC; now adding copying of box info to k from l
l = chem.amber.AmberParm('toluene_cyclohexane_10_500.prmtop', 'toluene_cyclohexane_10_500.inpcrd')
k = chem.load_file('toluene_cyclohexane_10_500_parmed.top')
k.box=l.box
k.positions = l.positions
sysl = l.createSystem( nonbondedMethod=app.PME, nonbondedCutoff=1.0 * u.nanometers, constraints=None )
sysk = k.createSystem( nonbondedMethod=app.PME, nonbondedCutoff=1.0 * u.nanometers, constraints=None )
conl = Context(sysl, VerletIntegrator(0.001))
conk = Context(sysk, VerletIntegrator(0.001))
conl.setPositions(l.positions)
conk.setPositions(l.positions)
print(chem.openmm.utils.energy_decomposition(l, conl))
print(chem.openmm.utils.energy_decomposition(k, conk))
