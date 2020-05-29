import numpy as np
import random
import FLScellclass as fls
import Th17CellClass1 as Th

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
            Th17_0 = int(input('Please input the number of Th17 cells that you would like to initialize:\n'))
            FLS_0 = int(input('Please input the number of FLS cells that you would like to initialize:\n'))

            # Piece to calculate the total number of cells that will be placed. might not need this, but leaving in for now
            cellCount = Th17_0 + FLS_0
            # %%
            # Here we prompt the user for the geometry of how they want the cells placed. We will likely mostly use either the None, Transwell, or Physiological
            distType = input(
                'Please specifiy how you would like the cells to be distributed in space. Options are: None (N), Random (R), Transwell (T), or Physiological (P), Co-culture (C):\n')

            if distType.lower() == 't':
                dim = '2d'  # dimension of locations

            elif distType.lower() == 'r':
                placement = random.shuffle(np.arange(l ** 3))
                FLSPlace = placement[:FLS_0]
                Th17Place = placement[FLS_0:Th17_0]
                dim = '3d'  # dimension of location matrix

            elif distType.lower() == 'n':
                Th17Place = np.zeros((Th17_0))
                FLSPlace = np.zeros((FLS_0))
                l = 1
                dim = '1d'

            elif distType.lower() == 'p':
                print('This spatial arrangement isnt supported yet')

            elif distType.lower() == 'c':
                print('This spatial arrangement isnt supported yet')

            else:
                print('This isnt a supported mode.')
            # %%
            if input('Would you like to load an initial set of cytokine profiles? Y/N:\n').lower() == 'y':
                il2 = np.zeros((l, l, l)) + float(input('Please specify the starting value of IL-2:\n'))
                il6 = np.zeros((l, l, l)) + float(input('Please specify the starting value of IL-6:\n'))
                gmcsf = np.zeros((l, l, l)) + float(input('Please specify the starting value of GM-CSF:\n'))
                il17 = np.zeros((l, l, l)) + float(input('Please specify the starting value of IL-17:\n'))
                il23 = np.zeros((l, l, l)) + float(input('Please specify the starting value of IL-23:\n'))
                il1b = np.zeros((l, l, l)) + float(input('Please specify the starting value of IL-1b:\n'))
            else:
                il2 = np.zeros((l, l, l))
                il6 = np.zeros((l, l, l))
                gmcsf = np.zeros((l, l, l))
                il17 = np.zeros((l, l, l))
                il23 = np.zeros((l, l, l))
                il1b = np.zeros((l, l, l))
            # %%
            tSteps = np.arange(int(input(
                'How long (in hours) would you like the simulation to run? Value sampling will occur every 12 hours and at the last hour:\n')) + 1)

    # %%
    print('Thank you! We are now completing initialization of the model. This may take a minute.\n')
    # This next bit is going to actually do the placing of the cells/initialization of the arrays that contain cells.
    if dim == '2d':
        Th17place = [1] * Th17_0 + [0] * ((l ** 2) - Th17_0)
        random.shuffle(Th17place)  # this creates a vector of size 1000 populated with zeroes and ones randomly
        FLSplace = [1] * FLS_0 + [0] * ((l ** 2) - FLS_0)
        random.shuffle(FLSplace)  # this creates a vector of size 1000 populated with zeroes and ones randomly
        Th17cellmat = np.zeros((l, l), dtype=Th.Th17cell) # initialize Th17 placement matrix
        FLScellmat = np.zeros((l, l), dtype=fls.FLScell) # initialize FLScell amtrix
        Th17cellplace = np.zeros((l, l))  # initialize Th17 reference matrix
        # initialize all cytokine matrices
        gmcsfmat = np.zeros((l, l, 2), dtype=float)
        il6mat = np.zeros((l, l, 2), dtype=float)
        il1bmat = np.zeros((l, l, 2), dtype=float)
        il17mat = np.zeros((l, l, 2), dtype=float)
        il6mat = np.zeros((l, l, 2), dtype=float)
        il8mat = np.zeros((l, l, 2), dtype=float)
        tnfamat = np.zeros((l, l, 2), dtype=float)


        for i in range(int(l ** (1 / 2))):
            Th17cellmat[i, :] = Th17place[int(i * l):int(i * l + (l))] # make matrix w zeroes and ones
            FLScellmat[i, :] = FLSplace[int(i * l):int(i * l + (l))] # " "
            Th17cellplace[i, :] = Th17place[int(i * l):int(i * l + (l))] # make matrix w zeroes and ones
            # for row in range():
            for element in range(len(Th17cellmat[i, :])):
                if Th17cellmat[i, element] == 1:
                    Th17cellmat[i, element] = Th.Th17cell(pos=[i, element, 0]) # populate matrix with cells at ones
                else:
                    Th17cellmat[i, element] = 0 # leave seroes as is.
            for element in range(len(FLScellmat[i, :])):
                if FLScellmat[i, element] == 1:
                    FLScellmat[i, element] = fls.FLScell( pos=[i, element, 1])  # populate matrix with cells at ones
                else:
                    FLScellmat[i, element] = 0  # leave seroes as is.
                            # this code has been tested and verified to create a matrix w random placement of th17s

            # print(Th17cellmat)

    elif dim == '1d':
        Th17place = [1] * Th17_0 + [0] * ((l) - Th17_0)
        random.shuffle(Th17place)  # this creates a vector of size 1000 populated with zeroes and ones randomly
        FLSplace = [1] * FLS_0 + [0] * ((l) - FLS_0)
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

    elif dim == '3d': # this is still in progress
        Th17place = [1] * Th17_0 + [0] * ((l ** 3) - Th17_0)
        random.shuffle(Th17place)  # this creates a vector of size 1000 populated with zeroes and ones randomly
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
    
    # This is the beginning of the secrete algorithm. Gets positions of all cells in Th17cellmat
    for i in range(len(Th17cellmat[1,:])):
        for j in range(len(Th17cellmat[:,i])):
            if Th17cellplace[i,j] == 1:
                loc = Th17cellmat[i,j].pos
        # This is the beginning of the secrete algorithm. Gets cells to sense and secrete cytokines
    for i in range(len(Th17cellmat[1,:])):
        for j in range(len(Th17cellmat[:,i])):
            if Th17cellplace[i,j] == 1:
                Thiscell = Th17cellmat[i,j]
                loc = Thiscell.pos
                loc = tuple(loc)
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
                il17mat[loc] = il17_secreted
                gmcsfmat[loc] = gmcsf_secreted
                # rnadom placement of proliferated cells algorithm
                if Thiscell.dblordie(il6_sensed, il23_sensed, il2_sensed, il7_sensed) == 'dbl':
                    # Defining cardinal direction matrix values
                    check_up = Th17cellplace[i,j+1]
                    check_down = Th17cellplace[i,j-1]
                    check_l = Th17cellplace[i-1,j]
                    check_r = Th17cellplace[i+1,j]
                    choose = random.randint(1, 4) # choosing random direction for proliferation
                    # UP   will need to improve this code so that it checks next direction after checking one direction
                    if choose == 1:
                        if check_up == 0:
                            Th17cellplace[i, j + 1] = 1
                            Th17cellmat[i, j + 1] = Th.Th17cell(pos = [i, j + 1, 0])
                    # DOWN
                    elif choose == 2:
                        if check_down == 0:
                            Th17cellplace[i, j - 1] = 1
                            Th17cellmat[i, j - 1] = Th.Th17cell(pos = [i, j - 1, 0])
                    # LEFT
                    elif choose == 3:
                        if check_l == 0:
                            Th17cellplace[i - 1, j] = 1
                            Th17cellmat[i - 1, j] = Th.Th17cell(pos = [i - 1, j, 0])
                    # RIGHT
                    elif choose == 4:
                        if check_r == 0:
                            Th17cellplace[i + 1, j] = 1
                            Th17cellmat[i + 1, j] = Th.Th17cell(pos = [i + 1, j, 0])
                    else:
                        x = 1  # This is where we will "push" cells to make space
