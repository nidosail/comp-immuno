import random
import numpy as np

# loc needs to be a tuple
def checkaround(cellplacemat, loc):
    dim = len(cellplacemat.shape)
    if dim == 1:
        x = loc[0]
        xpos = cellplacemat[x + 1]
        xneg = cellplacemat[x - 1]
    elif dim == 2:
        x = loc[0]
        y = loc[1]
        xpos = cellplacemat[x + 1, y]
        xneg = cellplacemat[x - 1, y]
        ypos = cellplacemat[x, y + 1]
        yneg = cellplacemat[x, y - 1]
    elif dim == 3:
        x = loc[0]
        y = loc[1]
        z = loc[2]
        xpos = cellplacemat[x + 1, y, z]
        xneg = cellplacemat[x - 1, y, z]
        ypos = cellplacemat[x, y + 1, z]
        yneg = cellplacemat[x, y - 1, z]
        zpos = cellplacemat[x, y, z + 1]
        zneg = cellplacemat[x, y, z - 1]
    else:
        return print('Dimension not supported')
    prolif_check = 'not prolif'  # for use in the while loop, to check if cell has proliferated into empty space
    # UP   will need to improve this code so that it checks next direction after checking one direction
    choose = np.arange(0, dim * 2)  # choosing random direction for proliferation
    random.shuffle(choose)
    for n in choose:
        if n == 0:
            if xpos == 0:
                if dim == 2:
                    prolifloc = [x + 1, y]

                if dim == 3:
                    prolifloc = [x + 1, y, z]

                if dim == 1:
                    prolifloc = [x + 1]

                prolif_check = 'prolif'
                break
        # DOWN
        elif n == 1:
            if xneg == 0:
                if dim == 2:
                    prolifloc = [x - 1, y]

                if dim == 3:
                    prolifloc = [x - 1, y, z]

                if dim == 1:
                    prolifloc = [x - 1]

                prolif_check = 'prolif'
                break
        # LEFT
        elif n == 2:
            if ypos == 0:
                if dim == 2:
                    prolifloc = [x, y + 1]

                if dim == 3:
                    prolifloc = [x, y + 1, z]

                prolif_check = 'prolif'
                break

        # RIGHT
        elif n == 3:
            if yneg == 0:
                if dim == 2:
                    prolifloc = [x, y - 1]

                if dim == 3:
                    prolifloc = [x, y - 1, z]

                prolif_check = 'prolif'
                break
        elif n == 4:
            if zpos == 0:

                if dim == 3:
                    prolifloc = [x, y, z + 1]

                prolif_check = 'prolif'
                break
        elif n == 5:
            if zneg == 0:

                if dim == 3:
                    prolifloc = [x, y, z - 1]

                prolif_check = 'prolif'
                break
    return [prolifloc, prolif_check]

def pushCell(cellplacemat, loc):
    dim = len(cellplacemat.shape)
    random.randint(0, dim * 2 - 1)
    matx = cellplacemat.shape[0]
    maty = cellplacemat.shape[1]
    if dim == 3:
        matz = cellplacemat.shape[2]
        zpos = cellplacemat[loc[0], loc[1], matz - loc[2]:]
        zneg = cellplacemat[loc[0], loc[1], abs(loc[2] - matz):]
    xpos = cellplacemat[matx - loc[0]:, loc[1], loc[2]]
    xneg = cellplacemat[abs(loc[0] - matx):, loc[1], loc[2]]
    ypos = cellplacemat[loc[0], maty - loc[1]:, loc[2]]
    yneg = cellplacemat[loc[0], abs(loc[1] - maty):, loc[2]]
    choose = np.arange(0,2*dim)

    def checkforpush(vec):
        pushnum = 0 # initialize variabe describing number of cells needed to push
        for n in vec:
            if n > 0:
                pushnum += 1 # for each cell in the row, add one to the number needed to be pushed
            else:
                break # when the first zero is encountered, break the loop
        return pushnum # returns the number of cells needed to be pushed
    for n in choose:
        if n == 0:
            if 0 in xpos:
                pushnum = checkforpush(xpos)
                direc = 'x positive'
        if n == 1:
            if 0 in xneg:
                pushnum = checkforpush(xpos)
                direc = 'x negative'
        if n == 2:
            if 0 in ypos:
                pushnum = checkforpush(xpos)
                direc = 'y positive'
        if n == 3:
            if 0 in yneg:
                pushnum = checkforpush(xpos)
                direc = 'y negative'
        if n == 4:
            if 0 in zpos:
                pushnum = checkforpush(xpos)
                direc = 'z positive'
        if n == 5:
            if 0 in zneg:
                pushnum = checkforpush(xpos)
                direc = 'z negative'
        else:
            print('The area has no more space for new cells')

    return [direc, pushnum] # returns direction of push for proliferation and how many cells needed to be pushed
