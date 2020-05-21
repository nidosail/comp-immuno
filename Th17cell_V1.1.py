class Th17cell:
    # some random sh!t that may come into use later, still tryna figure out the role of variables defined outside of
    # initialization method
    growth_factor = 1  # what is this?
    thresh2 = 70  # threshold of il2 required for cell proliferation
    thresh6 = 100  # threshold of il6 reuired for cell proliferation
    thresh7 = 100  # threshold of il23 required for cell prolif
    v23 = -3  # value that determines maximum enhancement to speed of prolif in the presence of  IL-23
    k23 = 100  # vale that determines sensitivite to IL-23
    size = 1  # number of voxels that this cell occupies

    # initialize variables, initializtion method, whatever u wanna call it. takes in position of form [x,y,z]
    # need to add numpy arrays for delayed cytokine secretion
    def __init__(self, pos=[0, 0, 0], k17=100, v17=200, kgm=100, vgm=200, delay=4, dblTmr=12, actTmr=24, dieTmr=36,
                 divNum=0):
        self.k17 = k17  # michaelis menten half concentration to max rate blah blah blah constant (units arbitrary rn)
        self.v17 = v17  # michaelis menten max rate. this is a made up value, here to hold up the skeleton of a model
        self.kgm = kgm
        self.vgm = vgm
        self.pos = pos  # position within the 3d box specified in cellmats 3rd,4th,and 5th arrays. eek!
        self.dil17 = [0] * delay  # initialize rate of il17 production for the current and next 3 timesteps. ie, dil17[0] would be the current rate of production, dil17[1] would be the rate of production one timestep in the future.
        self.dgmcsf = [0] * delay  # initialize rate of gmcsf production for the current and next 3 timesteps
        self.dblTmr = dblTmr
        self.actTmr = actTmr
        self.dieTmr = dieTmr
        self.divNum = divNum

    def secrete(self, il6, il1b):
        # cytokine rate according to michaelis menten kinetics with switching functions, needs work.
        # currently both cytokines would have the same rate, will look into this later.

        # we are going for a spooky scary skeleton of a model rn

        dil17_0 = self.dil17[0]  # This pulls the values of the rate of IL17 that will be secreted at the current timestep
        dgmcsf_0 = self.dgmcsf[0]  # This pulls the value of rate of GMCSF that will be secreted at the current timestep

        # this code updates the values stored for future secretion. I think we need to think more about this and maybe change the functions.
        if self.actTmr > 0:
            self.dil17 = self.dil17[1:] + [self.v17 * (il6 / (self.k17 + il6)) * self.v17 * (il1b / (self.k17 + il1b))]
            self.dgmcsf = self.dgmcsf[1:] + [
                self.vgm * (il6 / (self.kgm + il6)) * self.vgm * (il1b / (self.kgm + il1b))]
            self.actTmr = self.actTmr - 1
        else:
            self.dil17 = self.dil17[1:] + [0]
            self.dgmcsf = self.dil17[1:] + [0]
        # return the values of cytokines to be secreted
        return [dil17_0, dgmcsf_0]

    def dblOrDie(self, il6, il23=0, il2=0, il7=0):
        if self.dblTmr <= 0 and self.actTmr > 0 and self.divNum < 6 and (
                il6 >= self.thresh6 or il2 >= self.thresh2 or il7 >= self.thresh7):
            self.dblTmr = 12 + self.v23 * il23 / (il23 + self.k23)
            self.dieTmr = 36
            self.divNum = self.divNum + 1
            return 'dbl'
        elif self.dieTmr == 0:
            return 'die'
