import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import math

#Import Euler integrator for solving ODE system of chemical species inside the cells
from CellModeller.Integration.CLEulerIntegrator import CLEulerIntegrator

max_cells = 300000


def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, max_cells=max_cells, jitter_z=False)
    integ = CLEulerIntegrator(sim, 2, max_cells)

    # use this file for reg too
    regul = ModuleRegulator(sim)
    # Only biophys and regulation
    sim.init(biophys, regul, None, integ)

    # Specify the initial cell and its location in the simulation
    sim.addCell(cellType=0, pos=(0,0,0)) 

    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)

    sim.pickleSteps = 20



def init(cell):
    # Specify mean and distribution of initial cell size
    cell.targetVol = 3.0 + random.uniform(0.0,0.5)
    # Specify growth rate of cells
    
    cell.growthRate = 1.0
    
    # Specify initial concentration of chemical species
    cell.species[:] = [9,10]

def specRateCL():
    return '''
    const float a1 = 5.0f;
    const float a2 = 5.0f;
    
    float u = species[0];
    float v = species[1];
    
    rates[0] = a1 / (1 + v*v*v*v) - u;
    rates[1] = a2 / (1 + u*u*u*u) - v;
    '''

def update(cells):
    #Iterate through each cell and flag cells that reach target size for division
    for (id, cell) in cells.iteritems():
        #cell.color = [0.0, cell.species[5]/100, 0.0]
        cell.color = [cell.species[0]/5, cell.species[1]/5, 0]
        if cell.volume > cell.targetVol:
            a = 1
            cell.asymm = [a,1]
            cell.divideFlag = True

def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    d1.targetVol = 3.0 + random.uniform(0.0,0.5)
    d2.targetVol = 3.0 + random.uniform(0.0,0.5)

