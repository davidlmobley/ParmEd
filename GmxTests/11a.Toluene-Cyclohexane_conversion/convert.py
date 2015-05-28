#!/usr/bin/env python

import openmoltools.utils
import chemistry

in_prmtop = 'toluene_cyclohexane_10_500.prmtop'
in_crd = 'toluene_cyclohexane_10_500.inpcrd'

#Convert AMBER to GROMACS via acpype
openmoltools.utils.convert_via_acpype( 'mixture', in_prmtop, in_crd, out_top = 'toluene_cyclohexane_10_500_acpype.top', out_gro = 'toluene_cyclohexane_10_500_acpype.gro' )

#Convert AMBER to GROMACS via ParmEd
structure = chemistry.amber.AmberParm( in_prmtop, in_crd )
#Gromacs topology
gromacs_topology = chemistry.gromacs.GromacsTopologyFile.from_structure( structure )
#Write
gromacs_topology.write( 'toluene_cyclohexane_10_500_parmed.top' ) 
chemistry.gromacs.GromacsGroFile.write( gromacs_topology, 'toluene_cyclohexane_10_500_parmed.gro')

#Create a second AMBER coordinate file so it will have the same precision as the gro file
gromacs_topology2 = chemistry.load_file( 'toluene_cyclohexane_10_500_parmed.top' )
gro = chemistry.gromacs.GromacsGroFile.parse( 'toluene_cyclohexane_10_500_parmed.gro' )
gromacs_topology2.box = gro.box
gromacs_topology2.positions = gro.positions
amberparm = chemistry.amber.AmberParm.from_structure( gromacs_topology2 ) 
amberparm.write_parm( 'toluene_cyclohexane_10_500_parmed.prmtop' )
amberparm.write_rst7( 'toluene_cyclohexane_10_500_parmed.crd' )
