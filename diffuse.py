from itertools import product

def neighbors(index): # this is a function for finding the neighbor indices of a given index in a matrix for cytokine diffusion.
    N = len(index) # This has ben copy pasted lmao
    for relative_index in product((-1, 0, 1), repeat=N):
        if not all(i == 0 for i in relative_index):
            yield tuple(i + i_rel for i, i_rel in zip(index, relative_index))
