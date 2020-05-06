import numpy as np
import scipy as sp
import copy
import random

class FLScell:

    Thresh2 = 70  # threshold of il2 required for cell prolif
    Thresh6 = 100  # threshold of il6 reuired for cell proliferation
    Thresh7 = 100  # threshold of il23 required for cell prolif
    v23 = -3  # value that determines maximum enhancement to speed of prolif in the presence of  IL-23
    k23 = 100  # vale that determines sensitivite to IL-23
    size = 1  # number of voxels that this cell occupies

    # initialize variables, initializtion method, whatever u wanna call it. takes in position of form [x,y,z]
    def __init__(self, pos=[0, 0, 0], k8=100, v8=200, k6=100, v6=200, k1b = 100, v1b = 200, dil1b=np.zeros(4), dil6=np.zeros(4),
                 dil8=np.zeros(4), dblTmr=12, actTmr=0, dieTmr=36, divNum=0):
        # Michaelis-Menten constant initiation
        self.k8 = k8  # michaelis menten half concentration to max rate blah blah blah constant (units arbitrary rn)
        self.v8 = v8  # michaelis menten max rate. this is a made up value, here to hold up the skeleton of a model
        self.k6 = k6
        self.v6 = v6
        self.k1b = k1b
        self.v1b = v1b
        # other vars
        self.dil1b = dil1b  # initialize rate of il17 production for the current and next 3 timesteps. ie, dil17[0] would be the current rate of production, dil17[1] would be the rate of production one timestep in the future.
        self.dil6 = dil6  # initialize rate of gmcsf production for the current and next 3 timesteps
        self.dil8 = dil8
        self.pos = pos  # store the cells position. I think this wont be necessary.
        self.dblTmr = dblTmr
        self.actTmr = actTmr
        self.dieTmr = dieTmr
        self.divNum = divNum

    def secrete(self, gmcsf,il1b,tnfa):
        # cytokine rate according to michaelis menten kinetics with switching functions, needs work.
        # currently both cytokines would have the same rate, will look into this later.

        # we are going for a spooky scary skeleton of a model rn

        dil1b_0 = self.dil1b[0]  # This pulls the values of the rate of IL17 that will be secreted at the current timestep
        dil6_0 = self.dil6[0]  # This pulls the value of rate of GMCSF that will be secreted at the current timestep
        dil8_0 = self.dil8[0]

        # this code updates the values stored for future secretion. I think we need to think more about this and maybe change the functions.
        if self.actTmr > 0:
            self.dil1b = self.dil1b[1:] + [
                self.v1b * (gmcsf / (self.k1b + gmcsf)) * (il1b / (self.k1b + il1b)) * (tnfa / (self.k1b + tnfa))]
            self.dil6 = self.dil6[1:] + [
                self.v6 * (gmcsf / (self.k6 + gmcsf)) * (il1b / (self.k6 + il1b)) * (tnfa / (self.k6 + tnfa))]
            self.dil8 = self.dil8[1:] + [
                self.v8 * (gmcsf / (self.k8 + gmcsf)) * (il1b / (self.k8 + il1b)) * (tnfa / (self.k8 + tnfa))]
            self.actTmr = self.actTmr - 1
        else:
            self.dil17 = self.dil17[1:] + [0]
            self.dgmcsf = self.dil17[1:] + [0]
        # return the values of cytokines to be secreted
        return [dil1b_0, dil6_0, dil8_0]

    def dblOrDie(self, il6, il23=0, il2=0, il7=0):
        if self.dblTmr <= 0 and self.actTmr > 0 and divNum < 6 and (il6 >= self
        .Thresh6 or il2 >= self.Thresh2 or il7 >= self.Thresh7):
            self.dblTmr = 12 + self.v23 * il23 / (il23 + self.k23)
            self.dieTmr = 36
            self.divNum = self.divNum + 1
            return 'dbl'
        elif dieTmr == 0:
            return 'die'
