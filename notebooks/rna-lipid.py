import mdtraj as md
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from openmm import LangevinIntegrator
from openmm.app import PDBFile, ForceField, Simulation

# Load the RNA and Lipid Membrane System
pdb = PDBFile("rna_lipid_system.pdb")
forcefield = ForceField("amber99sbildn.xml", "lipid14.xml", "tip3p.xml")

# Create the system
topology = pdb.topology
system = forcefield.createSystem(topology, nonbondedMethod=app.PME,
                                 nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds)

# Define the integrator
integrator = LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 2.0*unit.femtoseconds)

# Set up the simulation
platform = mm.Platform.getPlatformByName("CUDA")  # Use GPU acceleration if available
properties = {"CudaPrecision": "mixed"}
simulation = Simulation(topology, system, integrator, platform, properties)
simulation.context.setPositions(pdb.positions)

# Minimize energy
print("Minimizing energy...")
simulation.minimizeEnergy()

# Run equilibration
equilibration_steps = 100000  # 100 ps
print("Running equilibration...")
simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
simulation.step(equilibration_steps)

# Run Production MD
production_steps = 500000  # 500 ps
print("Running production MD...")
simulation.reporters.append(app.DCDReporter("rna_lipid_md.dcd", 1000))
simulation.reporters.append(app.StateDataReporter("rna_lipid_md.log", 1000, step=True, potentialEnergy=True, temperature=True))
simulation.step(production_steps)

print("Simulation complete!")

