import random

def checkaround(cellplacemat, i, j):

    check_up = cellplacemat[i , j +1]
    check_down = cellplacemat[i , j -1]
    check_l = cellplacemat[ i -1 ,j]
    check_r = cellplacemat[ i +1 ,j]
    prolif_check = 'not proliferated'  # for use in the while loop, to check if cell has proliferated into empty space
    # UP   will need to improve this code so that it checks next direction after checking one direction
    while prolif_check != 'proliferated':
        choose = random.randint(1, 4)  # choosing random direction for proliferation
        # UP
        if choose == 1:
            if check_up == 0:
                cellplacemat[i, j + 1] = 1
                prolifloc = [i, j + 1]
                prolif_check = 'proliferated'
            else
                continue
        # DOWN
        elif choose == 2:
            if check_down == 0:
                cellplacemat[i, j - 1] = 1
                prolifloc = [i, j - 1]
                prolif_check = 'proliferated'
                break
            else
                continue
        # LEFT
        elif choose == 3:
            if check_l == 0:
                cellplacemat[i - 1, j] = 1
                prolifloc = [i - 1, j]
                prolif_check = 'proliferated'
                break
            else
                continue
        # RIGHT
        elif choose == 4:
            if check_r == 0:
                cellplacemat[i + 1, j] = 1
                prolifloc = [i + 1, j]
                prolif_check = 'proliferated'
                break
            else
                continue
    return prolifloc
