import random
import numpy as np
import vNNeighborCheck as vn
import Th17cellclass as Th
Th17_0 = 50
l = 10
Th17place = [1] * Th17_0 + [0] * ((l ** 2) - Th17_0)
random.shuffle(Th17place)  # this creates a vector of size 1000 populated with zeroes and ones randomly
Th17cellplace = np.zeros((l, l))
for i in range(int(l)):
    Th17cellplace[i, :] = Th17place[int(i * l):int(i * l + (l))] # make matrix w zeroes and ones

cells = np.where(Th17cellplace == 1) # finds the location of each integer 1 in the array
listofcells = list(zip(cells[0], cells[1])) # create a list out of the row ad column indices
print(neighbors(5,5))
loc = tuple(listofcells[25])

pushinfo = vn.pushCell(Th17cellplace, loc)
print(pushinfo)

mycell = Th.Th17cell(pos=[0, 0, 1])
cellpos = mycell.pos
mycell.pos = (2, 2, 2)
print(cellpos)
print(mycell.pos.__class__)
'''
for cell in listofcells:
    print(cell) # print individual indices
'''
'''
class Myclass:
    classnum = [2, 5, 9]
    def __init__(self):
        self.objnum = [0, 0, 0]
    def Mymethod(self,num):
        self.objnum = self.objnum + [num, num, num]



obj = Myclass()
obj.Mymethod(2)
print(obj.objnum)
'''
