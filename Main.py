import numpy as np
import random
import FLScell as fls
import Th17cellclass as Th
import vNNeighborCheck as vn
import diffuse as diff

def main():
    # First, greet the user. I want to make this a sort of user friendly interface to some degree,
    # so I'm starting out with that in mind. But, I am giving a hidden "skip" option in the beginning
    # for when we want to run a bunch of simulations, and it will basically skip straight to have you state how many
    # sims you want to run, and uploading a .txt file or something of the sort that contains your starting parameters

    if not input(
            'Hello, and welcome to the Agent-Based Cytokine-Driven model of Rheumatoid Arthritis, or "The ABCDs of RA". Press enter to begin.\n').lower() == 'skip':

        # This next section is again in place for user friendliness, but will mostly be in there for users who are want to be walked through how the model is initialized.
        if input(
                'Would you like to use the prompt interface to set initialization parameters, or would you prefer to use a .txt file to load parameters. For "prompt" enter "P" and for .txt enter "T"\n').lower() == 'p':
            # %%            #Load the default matrix size for now

            # Prompt asking if you want to change from the default lattice size
            if input(
                    'The default size is a square lattice of 100 voxel sides (1 mm^3). Would you like to change the default size? Y/N\n').lower() == 'y':
                l = int(input('Please put the integer number that represents a side of the lattice.\n'))
            else:
                l = 100

            # %%
            # Prompts to get the starting number of various cell types.
            Th17count = int(input('Please input the number of Th17 cells that you would like to initialize:\n'))
            FLScount = int(input('Please input the number of FLS cells that you would like to initialize:\n'))

            # Piece to calculate the total number of cells that will be placed. might not need this, but leaving in for now
            cellCount = Th17count + FLScount
            # %%
            # Here we prompt the user for the geometry of how they want the cells placed. We will likely mostly use either the None, Transwell, or Physiological
            distType = input(
                'Please specifiy how you would like the cells to be distributed in space. Options are: None (N), Random (R), Transwell (T), or Physiological (P), Co-culture (C):\n')

            if distType.lower() == 't':
                dim = 2  # dimension of locations

            elif distType.lower() == 'r':
                placement = random.shuffle(np.arange(l ** 3))
                FLSPlace = placement[:FLScount]
                Th17Place = placement[FLScount:Th17count]
                dim = 3  # dimension of location matrix

            elif distType.lower() == 'n':
                Th17Place = np.zeros((Th17count))
                FLSPlace = np.zeros((FLScount))
                l = 1
                dim = 1

            elif distType.lower() == 'p':
                print('This spatial arrangement isnt supported yet')

            elif distType.lower() == 'c':
                print('This spatial arrangement isnt supported yet')

            else:
                print('This isnt a supported mode.')
            # %%
            if input('Would you like to load an initial set of cytokine profiles? Y/N:\n').lower() == 'y':
                il2 = float(input('Please specify the starting value of IL-2:\n'))
                il6 = float(input('Please specify the starting value of IL-6:\n'))
                gmcsf = float(input('Please specify the starting value of GM-CSF:\n'))
                il17 = float(input('Please specify the starting value of IL-17:\n'))
                il23 = float(input('Please specify the starting value of IL-23:\n'))
                il1b = float(input('Please specify the starting value of IL-1b:\n'))
            else:
                il2 = 0
                il6 = 0
                gmcsf = 0
                il17 = 0
                il23 = 0
                il1b = 0
            # %%
            tSteps = np.arange(int(input(
                'How long (in hours) would you like the simulation to run? Value sampling will occur every 12 hours and at the last hour:\n')) + 1)

    # %%
    print('Thank you! We are now completing initialization of the model. This may take a minute.\n')
    # This next bit is going to actually do the placing of the cells/initialization of the arrays that contain cells.
    if dim == 2:
        Th17place = [1] * Th17count + [0] * ((l ** 2) - Th17count)
        random.shuffle(Th17place)  # this creates a vector of size 1000 populated with zeroes and ones randomly
        FLSplace = [1] * FLScount + [0] * ((l ** 2) - FLScount)
        random.shuffle(FLSplace)  # this creates a vector of size 1000 populated with zeroes and ones randomly
        Th17cellmat = np.zeros((l, l), dtype=Th.Th17cell) # initialize Th17 placement matrix
        FLScellmat = np.zeros((l, l), dtype=fls.FLScell) # initialize FLScell amtrix
        Th17cellplace = np.zeros((l, l))  # initialize Th17 reference matrix
        FLScellplace = np.zeros((l, l))  # initialize FLS reference matrix
        # initialize all cytokine matrices
        gmcsfmat = np.zeros((l, l, 2), dtype=float) + gmcsf
        il6mat = np.zeros((l, l, 2), dtype=float) + il6
        il1bmat = np.ones((l, l, 2), dtype=float) + il1b
        il17mat = np.zeros((l, l, 2), dtype=float) + il17
        il7mat = np.ones((l, l, 2), dtype=float)
        il8mat = np.zeros((l, l, 2), dtype=float)
        tnfamat = np.zeros((l, l, 2), dtype=float)
        il2mat = np.zeros((l, l, 2), dtype=float) + il2
        il23mat = np.zeros((l, l, 2), dtype=float) + il23


        for i in range(int(l)):
            Th17cellmat[i, :] = Th17place[int(i * l):int(i * l + (l))] # make matrix w zeroes and ones
            FLScellmat[i, :] = FLSplace[int(i * l):int(i * l + (l))] # " "
            Th17cellplace[i, :] = Th17place[int(i * l):int(i * l + (l))] # make matrix w zeroes and ones
            FLScellplace[i,:] = FLSplace[int(i * l):int(i * l + (l))]
            # for row in range():
            for element in range(len(Th17cellmat[i, :])):
                if Th17cellmat[i, element] == 1:
                    Th17cellmat[i, element] = Th.Th17cell(pos=[i, element, 0]) # populate matrix with cells at ones
                else:
                    Th17cellmat[i, element] = 0 # leave seroes as is.
            for element in range(len(FLScellmat[i, :])):
                if FLScellmat[i, element] == 1:
                    FLScellmat[i, element] = fls.FLScell(pos=[i, element, 1])  # populate matrix with cells at ones
                else:
                    FLScellmat[i, element] = 0  # leave seroes as is.
                            # this code has been tested and verified to create a matrix w random placement of th17s

            # print(Th17cellmat)

    elif dim == 1:
        Th17place = [1] * Th17count + [0] * ((l) - Th17count)
        random.shuffle(Th17place)  # this creates a vector of size 1000 populated with zeroes and ones randomly
        FLSplace = [1] * FLScount + [0] * ((l) - FLScount)
        random.shuffle(FLSplace)  # this creates a vector of size 1000 populated with zeroes and ones randomly
        Th17cellmat = np.zeros((l), dtype=Th.Th17cell)  # initialize Th17 placement matrix
        FLScellmat = np.zeros((l), dtype=fls.FLScell)  # initialize FLScell matrix

        for element in range(len(Th17cellmat)):
            if Th17cellmat[element] == 1:
                Th17cellmat[element] = Th.Th17cell(pos=[element, 0, 0])  # populate matrix with cells at ones
            else:
                Th17cellmat[element] = 0  # leave zeroes as is.
        for element in range(len(FLScellmat)):
            if FLScellmat[element] == 1:
                FLScellmat[element] = fls.FLScell(pos=[element, 0, 0])  # populate matrix with cells at ones
            else:
                FLScellmat[element] = 0  # leave zeroes as is.

    elif dim == 3: # this is still in progress
        Th17place = [1] * Th17count + [0] * ((l ** 3) - Th17count)
        random.shuffle(Th17place)  # this creates a vector of size 1000 populated with zeroes and ones randomly
        FLSplace = [1] * FLScount + [0] * ((l ** 3) - FLScount)
        Th17cellmat = np.zeros((l, l, l), dtype=Th.Th17cell) # initialize Th17 placement matrix
        FLScellmat = np.zeros((l, l, l), dtype=fls.FLScell) # initialize FLScell amtrix
        # in progress
        for i in range(l):
            for j in range(l):
                Th17cellmat[i, j, :] = Th17place[int(i * l):int(i * l + (l))] # make matrix w zeroes and ones
                FLScellmat[i, j, :] = FLSplace[int(i * l):int(i * l + (l))] # " "
                for element in range(l):
                    if Th17cellmat[i, element] == 1:
                        Th17cellmat[i, element] = Th.Th17cell(pos=[i, element, 0]) # populate matrix with cells at ones
                    else:
                        Th17cellmat[i, element] = 0 # leave seroes as is.
                for element in range(l):
                    if FLScellmat[i, element] == 1:
                        FLScellmat[i, element] = fls.FLScell( pos=[i, element, 0])  # populate matrix with cells at ones
                    else:
                        FLScellmat[i, element] = 0  # leave seroes as is.
                                # this code has been tested and verified to create a matrix w random placement of th17s
    for t in tSteps:
    # This is the beginning of the secrete algorithm. Gets cells to sense and secrete cytokines
    # Th17 cell dynamics
        for i in range(len(Th17cellmat[1,:])):
            for j in range(len(Th17cellmat[:,i])):
                if Th17cellplace[i,j] == 1:
                    Thiscell = Th17cellmat[i,j]
                    loc = tuple(Thiscell.pos)
                    # sensing cytokines
                    il6_sensed = il6mat[loc]
                    il1b_sensed = il1bmat[loc]
                    il7_sensed = il7mat[loc]
                    il2_sensed = il2mat[loc]
                    il23_sensed = il23mat[loc]
                    # secreting cytokines based on sensed cytokines
                    il17_gmcsf = Thiscell.secrete(il6_sensed, il1b_sensed)
                    # secreted cytokines
                    il17_secreted = il17_gmcsf[0]
                    gmcsf_secreted = il17_gmcsf[1]
                    # placing secreted cytokines in respective matrix
                    il17mat[loc] += il17_secreted
                    gmcsfmat[loc] += gmcsf_secreted
                    for place in list(diff.neighbors(loc)):
                        if max(place) <= max(gmcsfmat.shape) - 1 and place[2] <= gmcsfmat.shape[2] - 1:
                            il17mat[place] += il17_secreted
                            gmcsfmat[place] += gmcsf_secreted
                    # rnadom placement of proliferated cells algorithm
                    if Thiscell.dblordie(il6_sensed, il23_sensed, il2_sensed, il7_sensed) == 'dbl':
                        if 0 in Th17cellplace:
                            Th17count += 1
                        print(Th17count)
                        # This is what happens when cells proliferate
                        # Defining cardinal direction matrix values
                        prolifinf = vn.checkaround(Th17cellplace, i, j) # this uses the written function to check all four directions in a separate file]
                        if prolifinf[1] == 'prolif':
                            prolifloc = prolifinf[0]
                            # prolifloc needs to be tuple to be referenced as the indices of a matrix. places Th17cell there
                            # need to make improvements to prolifcheck file
                            Th17cellmat[tuple(prolifloc)] = Th.Th17cell(pos=prolifloc)
                            Th17cellplace[tuple(prolifloc)] = 1
                        else:
                            pushloc = vn.pushCell(Th17cellplace) # This returns two arrays: a "unit vector" of size 1x3 that defines direction
                            # and another vector of size 1x3 which defines the length to be pushed, to define the vector of pushing
                            push = tuple(np.add(loc,pushloc[1])) # a variable which defines the end location of the array of pushed cells
                            newloc = tuple(np.add(loc, pushloc[0])) # this variable defines the location of the new cell
                            for cell in Th17cellmat[newloc[0]:push[0],newloc[1]:push[1],newloc[2]:push[2]]:
                                cellpos = tuple(cell.pos)
                                cellnextpos = tuple(np.add(cellpos,pushloc[0]))
                                cell.pos = cellnextpos
                                Th17cellmat[cellnextpos] = Th17cellmat[cellpos]
                                Th17cellplace[cellnextpos] = 1
                            Th17cellmat[newloc] = Th.Th17cell(pos=newloc)

                    elif Thiscell.dblordie(il6_sensed, il23_sensed, il2_sensed, il7_sensed) == 'die':
                        # This is what happens when cells die
                        Th17cellplace[loc] = 0
                        Th17cellmat[loc] = 0
        # FLS cell dynamics
        for i in range(len(FLScellmat[1,:])):
            for j in range(len(FLScellmat[:,i])):
                if FLScellplace[i,j] == 1:
                    Thisflscell = FLScellmat[i,j]
                    loc = tuple(Thisflscell.pos)
                    # sensing cytokines
                    il6_sensed = il6mat[loc]
                    il1b_sensed = il1bmat[loc]
                    il7_sensed = il7mat[loc]
                    il2_sensed = il2mat[loc]
                    il23_sensed = il23mat[loc]
                    il17_sensed = il17mat[loc]
                    gmcsf_sensed = gmcsfmat[loc]
                    tnfa_sensed = tnfamat[loc]
                    # secreting cytokines based on sensed cytokines
                    il1b_6_8 = Thisflscell.secrete(gmcsf_sensed, il1b_sensed, tnfa_sensed)
                    # secreted cytokines
                    il1b_secreted = il1b_6_8[0]
                    il6_secreted = il1b_6_8[1]
                    il8_secreted = il1b_6_8[2]
                    # placing secreted cytokines in respective matrix
                    il1bmat[loc] += il1b_secreted  # secreting cytokines at the cells locations.
                    il6mat[loc] += il6_secreted
                    il8mat[loc] += il8_secreted
                    for place in list(diff.neighbors(loc)): # This calculates the diffusion of the cytokines and places them
                        if max(place) <= max(gmcsfmat.shape) - 1 and place[2] <= gmcsfmat.shape[2] - 1:
                            il1bmat[place] += il1b_secreted # at the diffusion indices
                            il6mat[place] += il6_secreted
                            il8mat[place] += il8_secreted
                    # rnadom placement of proliferated cells algorithm
                    if Thisflscell.dblordie(il6_sensed, il23_sensed, il2_sensed, il7_sensed) == 'dbl':
                        if 0 in FLScellplace:
                            FLScount += 1
                        print(FLScount)
                        # This is what happens when cells proliferate
                        # Defining cardinal direction matrix values
                        prolifinf = vn.checkaround(FLScellplace, loc, dim) # this uses the written function to check all four directions in a separate file]
                        if prolifinf[1] == 'prolif':
                            prolifloc = prolifinf[0]
                            # prolifloc needs to be tuple to be referenced as the indices of a matrix. places Th17cell there
                            # need to make improvements to prolifcheck file
                            FLScellmat[tuple(prolifloc)] = fls.FLScell(pos=prolifloc)
                            FLScellplace[tuple(prolifloc)] = 1
                        else:
                            pushinfo = vn.pushCell(FLScellplace, loc)

                    elif Thisflscell.dblordie(il6_sensed, il23_sensed, il2_sensed, il7_sensed) == 'die':
                        # This is what happens when cells die
                        FLScellplace[loc] = 0
                        FLScellmat[loc] = 0
