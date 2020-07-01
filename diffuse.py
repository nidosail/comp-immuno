import numpy as np

def around6(loc):
    listoflocs = np.array([loc])
    for num in loc:
        np.append(listoflocs,[[num]])
