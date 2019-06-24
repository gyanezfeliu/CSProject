import matplotlib
import sys
import os
import math
import numpy as np
import pandas as pd
sys.path.append('.')
import cPickle
import CellModeller
import json
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# returns an array that contains IDs for cells that are in the area (ymax-ymin)*diameter
def cells_in_center_section(cs, it, ymin, ymax):
    pos_x = []
    pos_y = []

    for it in cs:
        x = cs[it].pos[0]
        y = cs[it].pos[1]
        pos_x.append(x)
        pos_y.append(y)

    # min and max of x and y
    xmin = math.ceil(min(pos_x))
    xmax = math.floor(max(pos_x))

    # selects cells that are in x and y range
    cells_range = []
    for c in cs:
        if cs[c].pos[0] >= xmin and cs[c].pos[0] <= xmax:
            if cs[c].pos[1] >= ymin and cs[c].pos[1] <= ymax:
                cells_range.append(c)
    cells_range = np.array(cells_range)
    
    return cells_range

# returns cell states dictionary from a pickle
def get_cs(path, file_name):
    data = cPickle.load(open(path+file_name,'r'))
    cs= data['cellStates']
    it = iter(cs)
    return cs, it

# takes two consecutive pickles and returns an array that contains the IDs of the cells that
# are in both pickles
def get_matchcells(cells_range1, cells_range2):
    match_1to2 = np.zeros(len(cells_range1))
    for i, cell in enumerate(cells_range1):
        if cell in cells_range2:
            match_1to2[i] = cell
        else:
            match_1to2[i] = 0   
    match_index = cells_range1[np.where(match_1to2 != 0)]
    return match_index

# computes effective growth rate for cells in "match_index" between two time steps
# stored in cs1 and cs2
def calc_effGrowthRates(match_index, cs1, cs2):
    x = []
    effGrowthRate = []
    for index in match_index:
        # considers just colony's right side
        if cs2[index].pos[0] >= 0:
            effGrowthRate.append((cs2[index].length - cs1[index].length) / cs1[index].length)
            x.append(cs2[index].pos[0])
    effGrowthRate = np.array(effGrowthRate)
    return x, effGrowthRate

# it transforms a continuous domain space into a discrete domain space
def discretise_domain(x, x_gt0=False):
    x_range = np.zeros(len(x))
    for i, val in enumerate(x):
        x_range[i] = math.floor(val)
    x_prom = np.unique(x_range)
    # if True considers just colony's right side
    if x_gt0:
        x_prom = x_prom[np.where(x_prom >=0)]
    return x_range, x_prom

# computes average effective growth rate for a discretized domain
def avg_measure(x_prom, x_range, measure):
    measure_prom = np.zeros(len(x_prom))
    for index, x in enumerate(x_prom):
        ind = np.where(x_range == x)[0]
        ## Average growth rate growth rate for position i in domain X
        measure_prom[index] = np.mean(measure[ind])
    
    return measure_prom

# computes drift velocity for cells that are present in two consecutive pickles
def calc_driftVelocity(match_index, cs1, cs2):
    x = [] 
    driftVelocity = []
    for ID in match_index:
        driftVelocity.append(np.sqrt(cs2[ID].pos[0]**2 + cs2[ID].pos[1]**2) - np.sqrt(cs1[ID].pos[0]**2 + cs1[ID].pos[1]**2))
        x.append(cs2[ID].pos[0])
    driftVelocity = np.array(driftVelocity)
    return x, driftVelocity

# returns effective growth rate
def get_effGrowthRates(path, files, ymin=-10, ymax=10):
    x_prom = []
    effgrs = []
    for i in range(len(files)-2):
        print(i)
        ### Generate IDs of the cells present in steps t and t+1
        cs2, it2 = get_cs(path, files[i+1])
        cs1, it1 = get_cs(path, files[i])
        cells_range1 = cells_in_center_section(cs1, it1, ymin, ymax)
        cells_range2 = cells_in_center_section(cs2, it2, ymin, ymax)
        ### Get cells that are present in step t+1 
        match_index = get_matchcells(cells_range1, cells_range2)
        ### Calculate effGrowthRate
        x, effGrowthRate = calc_effGrowthRates(match_index, cs1, cs2)
        ### Discretise radius domain
        x_range, x_prom = discretise_domain(x, x_gt0=False)        
        ### Extracs data from cells thar are in range
        effgrs = avg_measure(x_prom, x_range, effGrowthRate)
    return x_prom, effgrs

# returns drift velocity
def get_DriftVelocity(path, files, ymin=-10, ymax=10):
    x_prom = []
    drftvel = []
    for i in range(len(files)-2):
        print(i)
        ### Generate IDs of the cells present in steps t and t+1
        cs2, it2 = get_cs(path, files[i+1])
        cs1, it1 = get_cs(path, files[i])
        cells_range1 = cells_in_center_section(cs1, it1, ymin, ymax)
        cells_range2 = cells_in_center_section(cs2, it2, ymin, ymax)
        ### Get cells that are present in step t+1 
        match_index = get_matchcells(cells_range1, cells_range2)
        ### Calculate driftVelocity        
        x, driftVelocity = calc_driftVelocity(match_index, cs1, cs2)
        ### Discretise radius domain
        x_range, x_prom = discretise_domain(x, x_gt0=True)        
        ### Extracs data from cells thar are in range
        drftvel = avg_measure(x_prom, x_range, driftVelocity)    
    return x_prom, drftvel

# tracks single cell and stores position, color, current daughter and step number
def track_colorPos(path, files):
    data = cPickle.load(open(path+files[len(files)-1],'r'))
    ln = data['lineage']
    dot = np.array(ln.keys())
    par = np.array(ln.values())

    # states
    pos = []
    color = []
    daugh= []
    stepNum = []

    # First daughter
    daughter = 2
    for i, file_name in enumerate(files):
        print(i)
        data = cPickle.load(open(path+file_name, 'r'))
        cs = data['cellStates']
        cs_keys = np.array(cs.keys())
        # if daughter is not in picke we have to find the first daughter of current daughter
        if daughter not in cs_keys:
            while daughter not in cs_keys:
                daughter = dot[np.where(par==daughter)[0][0]]
        # extract states
        pos.append(cs[daughter].pos)
        color.append(cs[daughter].color)
        daugh.append(daughter)
        stepNum.append(data['stepNum'])
    df = pd.DataFrame(data={'position': pos, 'color': color, 'daughter': daugh, 'stepNum': stepNum})
    return df

def get_fluorescence_levels(path, files):
    database = []
    for f in files:    
        cs, it = get_cs(path, f)
        cells_range = cells_in_center_section(cs, it, -10, 10)

        # extracts coordinates and fluorescence levels that are in the range
        x_range = np.zeros(len(cells_range))
        y_range = np.zeros(len(cells_range))
        R_range = np.zeros(len(cells_range))
        G_range = np.zeros(len(cells_range))
        B_range = np.zeros(len(cells_range))

        for i, c in enumerate(cells_range):
            x_range[i] = math.floor(cs[c].pos[0])
            y_range[i] = cs[c].pos[1]
            R_range[i] = cs[c].color[0]
            G_range[i] = cs[c].color[1]    
            B_range[i] = cs[c].color[2]

        # for unique values of x sums all R, G y B that are in the range and gets the average.
        x_prom = np.unique(x_range)
        R_prom = np.zeros(len(x_prom))
        G_prom = np.zeros(len(x_prom))
        B_prom = np.zeros(len(x_prom))
        for index, x in enumerate(x_prom):
            ind = np.where(x_range == x)
            ## Average fluorescence R, G, B for position i in domain X
            R_prom[index] = np.mean(R_range[ind])
            G_prom[index] = np.mean(G_range[ind])
            B_prom[index] = np.mean(B_range[ind])

        database.append([int(f[5:10]), list(x_prom), list(R_prom), list(G_prom), list(B_prom)])
        
    return database
    
def wave_time(T, lam, r0, ri, t):
    return (np.sin( ((2*np.pi) / T) * (t - ((r0 - ri) / lam) * T )))**2

def ri_update(v0, ri, dt, lam):
    vi = vi_calc(ri, lam, v0, lam_scale)
    return (ri + dt * vi)

def vi_calc(ri, lam, v0, lam_scale):
    return v0
    #return (np.tanh( (ri / (lam_scale*lam))) + 1 ) * v0

def make_circles(ri, phi_list):
    circle_x = [ri * np.sin(phi) for phi in phi_list]
    circle_y = [ri * np.cos(phi) for phi in phi_list]

    return circle_x, circle_y