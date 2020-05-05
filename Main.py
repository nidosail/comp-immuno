# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 20:09:48 2020

@author: davem
"""
# Lets make a function to run the simulation

import numpy as np
import random

# timestep vector (for experimental purposes for now) 48 hours? minutes? seconds? tacos? im hungry for tacos
def main():
    if not input('Hello, and welcome to the Agent Based Cytokine Driven model of Rheumatoid Arthritis, or "ABCDs of RA". Press enter to begin.\n').lower() == "skip":
        if input('Would you like to use the prompt interface to set initialization parameters, or would you prefer to use a .txt file to load parameters. For "prompt" enter "P" and for .txt enter "T"\n').lower() == 'p':
        
            lsize = 100
            if input('The default size is a square lattice of 100 voxel sides (1 mm^3). Would you like to change the default size? Y/N\n').lower() == 'y':
                lsize = int(input('Please put the integer number that represents a side of the lattice.\n'))
            
            Th17_0 = int(input('Please input the number of Th17 cells that you would like to initialize:\n'))
            FLS_0 = int(input('Please input the number of FLS cells that you would like to initialize:\n'))
            
            cellCount = Th17_0 + FLS_0
            
            distType = input('Please specifiy how you would like the cells to be distributed in space. Options are: None (N), Random (R), Planar (P).:\n')
            
            if distType.lower() == 'r':
                placement = random.shuffle(np.arange(lsize**3))[:cellCount]
            
            
    

