# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 20:08:16 2020

@author: davem
"""


import numpy as np
import scipy as sp
import copy
import random

# begin class code
class Th17cell:
    # some random sh!t that may come into use later, still tryna figure out the role of variables defined outside of
    # initialization method
    growth_factor = 1 #what is this?
    thresh2 = 70 #threshold of il2 required for cell proliferation
    thresh6 = 100 #threshold of il6 reuired for cell proliferation
    thresh7 = 100 #threshold of il23 required for cell prolif
    v23 = -3 #value that determines maximum enhancement to speed of prolif in the presence of  IL-23
    k23 = 100 #vale that determines sensitivite to IL-23
    size = 1 # number of voxels that this cell occupies
    # initialize variables, initializtion method, whatever u wanna call it. takes in position of form [x,y,z]
    def __init__(self, pos = [0,0,0], k17 = 100, v17=200, kgm = 100, vgm = 200, delay = 4, dblTmr = 12, actTmr = 24, dieTmr = 36, divNum = 0):
        self.k17 = k17 # michaelis menten half concentration to max rate blah blah blah constant (units arbitrary rn)
        self.v17 = v17 # michaelis menten max rate. this is a made up value, here to hold up the skeleton of a model
        self.kgm = kgm
        self.vgm = vgm
        self.pos = pos # position within the 3d box specified in cellmats 3rd,4th,and 5th arrays. eek!
        self.dil17 = [0]*delay # initialize rate of il17 production for the current and next 3 timesteps. ie, dil17[0] would be the current rate of production, dil17[1] would be the rate of production one timestep in the future.
        self.dgmcsf = [0]*delay # initialize rate of gmcsf production for the current and next 3 timesteps 
        self.pos = pos #store the cells position. I think this wont be necessary.
        self.dblTmr = dblTmr
        self.actTmr = actTmr
        self.dieTmr = dieTmr
        self.divNum = divNum
        
    def secrete(self, il6, il1b):
        # cytokine rate according to michaelis menten kinetics with switching functions, needs work.
        # currently both cytokines would have the same rate, will look into this later.

        # we are going for a spooky scary skeleton of a model rn
        
        dil17_0 = self.dil17[0] #This pulls the values of the rate of IL17 that will be secreted at the current timestep
        dgmcsf_0 = self.dgmcsf[0] #This pulls the value of rate of GMCSF that will be secreted at the current timestep
        
        #this code updates the values stored for future secretion. I think we need to think more about this and maybe change the functions.
        if self.actTmr > 0:
            self.dil17 = self.dil17[1:] + [self.v17*(il6/(self.k17 + il6)) * self.v17*(il1b/(self.k17 + il1b))]
            self.dgmcsf = self.dgmcsf[1:] + [self.vgm*(il6/(self.kgm + il6)) * self.vgm*(il1b/(self.kgm + il1b))]
            self.actTmr = self.actTmr - 1
        else:
            self.dil17 = self.dil17[1:] + [0]
            self.dgmcsf = self.dil17[1:] + [0]
        # return the values of cytokines to be secreted
        return [dil17_0, dgmcsf_0]
    
    def dblOrDie(self, il6, il23 = 0, il2 = 0, il7 = 0):
        if self.dblTmr <= 0 and self.actTmr > 0 and self.divNum < 6 and (il6 >= self.thresh6 or il2 >= self.thresh2 or il7 >= self.thresh7):
            self.dblTmr = 12 + self.v23*il23/(il23 + self.k23)
            self.dieTmr = 36
            self.divNum = self.divNum + 1
            return 'dbl'
        elif self.dieTmr == 0:
            return 'die'
#%%
#Function designed to pick a position in a Moores neighborhood about a certain point. 
def pos_picker(pos, dim = '2d', cellMat):
    

#%%
#This will end up being the function that runs the actual simulation and kinda just wraps everything together.
def main():
    #First, greet the user. I want to make this a sort of user friendly interface to some degree, 
    #so I'm starting out with that in mind. But, I am giving a hidden "skip" option in the beginning 
    #for when we want to run a bunch of simulations, and it will basically skip straight to have you state how many
    #sims you want to run, and uploading a .txt file or something of the sort that contains your starting parameters

    if not input('Hello, and welcome to the Agent-Based Cytokine-Driven model of Rheumatoid Arthritis, or "The ABCDs of RA". Press enter to begin.\n').lower() == 'skip':
        
        #This next section is again in place for user friendliness, but will mostly be in there for users who are want to be walked through how the model is initialized.
        if input('Would you like to use the prompt interface to set initialization parameters, or would you prefer to use a .txt file to load parameters. For "prompt" enter "P" and for .txt enter "T"\n').lower() == 'p':
#%%            #Load the default matrix size for now
            
            #Prompt asking if you want to change from the default lattice size
            if input('The default size is a square lattice of 100 voxel sides (1 mm^3). Would you like to change the default size? Y/N\n').lower() == 'y':
                l = int(input('Please put the integer number that represents a side of the lattice.\n'))
            else:
                l = 100
 
#%%
           #Prompts to get the starting number of various cell types.
            Th17_0 = int(input('Please input the number of Th17 cells that you would like to initialize:\n'))
            FLS_0 = int(input('Please input the number of FLS cells that you would like to initialize:\n'))
            
            #Piece to calculate the total number of cells that will be placed. might not need this, but leaving in for now
            cellCount = Th17_0 + FLS_0
#%%            
            #Here we prompt the user for the geometry of how they want the cells placed. We will likely mostly use either the None, Transwell, or Physiological
            distType = input('Please specifiy how you would like the cells to be distributed in space. Options are: None (N), Random (R), Transwell (T), or Physiological (P), Co-culture (C):\n')
                
            if distType.lower() == 't':
                loc_dim = '2d' #dimension of locations
                
            elif distType.lower() == 'r':
                placement = random.shuffle(np.arange(l**3))
                FLSPlace = placement[:FLS_0]
                Th17Place = placement[FLS_0:Th17_0]
                dim = '3d' #dimension of location matrix
                
            elif distType.lower() == 'n':
                Th17Place = np.zeros((Th17_0))
                FLSPlace = np.zeros((FLS_0))
                l = 1
                dim = '0d'
                
            elif distType.lower() == 'p':
                print('This spatial arrangement isnt supported yet')
                
            elif distType.lower() == 'c':
                print('This spatial arrangement isnt supported yet')
                
            else:
                print('This isnt a supported mode.')
#%%            
            if input('Would you like to load an initial set of cytokine profiles? Y/N:\n').lower() == 'y':
                il2 = np.zeros((l,l,l)) + float(input('Please specify the starting value of IL-2:\n'))
                il6 = np.zeros((l,l,l)) + float(input('Please specify the starting value of IL-6:\n'))
                gmcsf = np.zeros((l,l,l)) + float(input('Please specify the starting value of GM-CSF:\n'))
                il17 = np.zeros((l,l,l)) + float(input('Please specify the starting value of IL-17:\n'))
                il23 = np.zeros((l,l,l)) + float(input('Please specify the starting value of IL-23:\n'))
                il1b = np.zeros((l,l,l)) + float(input('Please specify the starting value of IL-1b:\n'))
            else:
                il2 = np.zeros((l,l,l))
                il6 = np.zeros((l,l,l))
                gmcsf = np.zeros((l,l,l))
                il17 = np.zeros((l,l,l))
                il23 = np.zeros((l,l,l))
                il1b = np.zeros((l,l,l))
#%%                
            tSteps = np.arange(int(input('How long (in hours) would you like the simulation to run? Value sampling will occur every 12 hours and at the last hour:\n')) + 1)

#%%           
    print('Thank you! We are now completing initialization of the model. This may take a minute.\n')
    #This next bit is going to actually do the placing of the cells/initialization of the arrays that contain cells. 
    if dim == '2d':
        Th17place = [1] * Th17_0 + [0] * ((l ** 2) - Th17_0)
        random.shuffle(Th17place)  # this creates a vector of size 10000, or l^2 populated with zeroes and ones randomly
        FLSplace = [1] * Th17_0 + [0] * ((l ** 2) - FLS_0)
        random.shuffle(FLSplace)  # this creates a vector of size 10000, or l^2 populated with zeroes and ones randomly
        Th17cellmat = np.zeros((l, l), dtype=Th17cell) # initialize Th17 placement matrix
        FLScellmat = np.zeros((l, l), dtype=Th17cell) # initialize FLScell amtrix

        for i in range(int(l)):
            Th17cellmat[i, :] = Th17place[int(i * l):int(i * l + (l))] # make matrix w zeroes and ones
            FLScellmat[i, :] = FLSplace[int(i * l):int(i * l + (l))] # " "
            # for row in range():
            for element in range(len(Th17cellmat[i, :])):
                if Th17cellmat[i, element] == 1:
                    Th17cellmat[i, element] = Th17cell(pos=[i, element, 0]) # populate matrix with cells at ones
                else:
                    Th17cellmat[i, element] = 0 # leave seroes as is.
            for element in range(len(FLScellmat[i, :])):
                if FLScellmat[i, element] == 1:
                    FLScellmat[i, element] = FLScell( pos=[i, element, 0])  # populate matrix with cells at ones
                else:
                    FLScellmat[i, element] = 0  # leave seroes as is.
                            # this code has been tested and verified to create a matrix w random placement of th17s
            
    elif dim == '3d':
        Th17s = np.zeros((l,l,l),dtype = Th17cell)
        FLSs = np.zeros((l,l,l),dtype = FLScell)
        for i in Th17Place:
            zpos = i//l**2
            ypos = (i - zpos*l**2)//l
            xpos = (i - zpos*l**2 - ypos*l)
            Th17s[xpos,ypos,zpos] = Th17cell(pos = [xpos,ypos,zpos])
        for i in FLSPlace:
            zpos = i//l**2
            ypos = (i - zpos*l**2)//l
            xpos = (i - zpos*l**2 - ypos*l)
            FLSs[xpos,ypos,zpos] = FLScell(pos = [xpos,ypos,zpos])
        
    cellMat = np.zeros((l,l,l))
        
        
#%%            
    print('Beginning simulation\n')
    
    dgmcsf = np.zeros((l,l,l))
    dil6 = np.zeros((l,l,l))
    dil17 = np.zeros((l,l,l))
    dil1b = np.zeros((l,l,l))
    
    for t in tSteps:
        print('On timestep:' + str(t) + '\n')
        #This loop will determine rates of cytokine production for solving the PDE
        for cell in Th17s:
            
            if not cell == 0:
                prolif = cell.dblOrDie(il6 = il6, il23 = il23, il2 = il2)
                cellMat(Th17 = 0)
                if prolif == 'die':
                    Th17s[cell.pos[0], cell.pos[1], cell.pos[2]] = 0
                    
                elif prolif == 'dbl':
                    [dil17[cell.pos[0], cell.pos[1], cell.pos[2]], dgmcsf[cell.pos[0], cell.pos[1], cell.pos[2]]] = cell.secrete(il6, il1b)
                    divLocs = shuffle(np.arange
                else:
                    [dil17[cell.pos[0], cell.pos[1], cell.pos[2]], dgmcsf[cell.pos[0], cell.pos[1], cell.pos[2]]] = cell.secrete(il6, il1b)
                
                
            
            
                
        for cell in FLSs:
            if not cell == 0:
                
                [dil6[cell.pos[0], cell.pos[1], cell.pos[2]]] = cell.secrete(gmcsf, il1b)
        
        
                
        
            
            
            
                
                
                
#%%
#csmajors#arechumps
#pythonfun
#chemewho?
#ihopedaveapproves
#daveisproud
#hashtag
