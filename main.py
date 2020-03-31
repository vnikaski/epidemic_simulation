import numpy as np
from random import randint, random
from matplotlib import pylab
from matplotlib.animation import FuncAnimation
import argparse


def update_neighbours_for_cell(map: np.array, direction: str, i: int, j: int, r: int):
    """
    Updates number of Moore's neighbours in the distance of r from the cell map[i,j]
    :param map: map of states
    :param direction: 'up', 'down', 'right', 'left'
    :param i: row of the cell
    :param j: column of the cell
    :param r: radius of Moore's neighbourhood
    :return: updated map: np.array
    """

    a = 0 #sum of infected neighbours in given direction
    for k in range(r):
        b = k #parameter needed in the while loop to check for the edges of the map
        c = k #same as above

        if direction == 'up':
            while j-b < 0:
                b -= 1
            while j+c+2 > len(map):
                c -= 1
            a = sum(map[i,j-b:j+c+2,0]==1)
        elif direction == 'down':
            while j-b-1 < 0:
                b -= 1
            while j+c+1 > len(map):
                c -= 1
            a = sum(map[i, j-b-1:j+c+1, 0]==1)
        elif direction == 'left':
            while i - b - 1 < 0:
                b -= 1
            while i + c + 1 > len(map):
                c -= 1
            a = sum(map[i-b-1:i+c+1, j, 0]==1)
        elif direction == 'right':
            while i-b < 0:
                b -= 1
            while i+c+2 > len(map):
                c -= 1
            a = sum(map[i-b:i+c+2, j, 0]==1)

        map[i,j,1] += a
    return map

def update_neighbours(map: np.array, r: int):
    """
    Goes through all of the map to update neighbours in every direction
    :param map: np.array map of states
    :param r: radius of infection
    :return: updated map np.array
    """
    for i in range(len(map)):
        for j in range(len(map)):
            map = update_neighbours_for_cell(map, 'up', i, j, r)
            map = update_neighbours_for_cell(map, 'right', i, j, r)
            map = update_neighbours_for_cell(map, 'down', i, j, r)
            map = update_neighbours_for_cell(map, 'left', i, j, r)
    return map

def main(N: int, k: int, p_w: float, p_z: float, M: int, r: int = 1):
    """
    Creates simulation of a spreading infection on a square map. Each cell is in one of the three states:
    0 - healthy, capable of getting infected
    1 - infected, can spread the infection
    2 - cured, no longer spreading, can't get infected

    :param N: size of the edge of  the square
    :param k: number of first randomly infected cells
    :param p_w: probability of curing the infection by an infected cell per epoch
    :param p_z: probability of getting the infection an infected neighbour cell (changes with the number of infected neighbours)
    :param M: number of epochs
    :param r: radius of spreadage
    """

    map = np.zeros((N,N,2)) #creating map; every cell has two dimensions: [state, number_of_infected_neighbours]
    while k > 0: #choosing randomly k infected people
        i = randint(0, N-1)
        j = randint(0, N-1)
        if map[i,j,0] == 0:
            map[i,j,0] = 1
            k -= 1
    map = update_neighbours(map, r) #updating infecting neighbours after random infection
    count = {0: [sum(sum(map[:, :, 0] == 0))], 1: [sum(sum(map[:, :, 0] == 1))], 2: [sum(sum(map[:, :, 0] == 2))]}


    #preparing for data storage needed for the animation
    maps = np.zeros((N, N, M))
    maps[:, :, 0] = map[:, :, 0]


    for e in range(M): #iterating through epochs
        for i in range(N): #going through rows of the map; i = row in
            for j in range(N):#going through columns of the map; j = column in
                if map[i,j,0] == 0 and map[i,j,1]>0 and random() < 1-(1-p_z)**map[i,j,1]: #trying to infect cell with probability = 1-(1-p_z)
                    map[i,j,0] = 1
                elif map[i,j,0] == 1 and random() < p_w: #trying to heal infected cell
                    map[i,j,0] = 2
        update_neighbours(map, r)

        #counting epoch stats
        count[0].append(sum(sum(map[:, :, 0] == 0)))
        count[1].append(sum(sum(map[:, :, 0] == 1)))
        count[2].append(sum(sum(map[:, :, 0] == 2)))

        #drawing and saving heatmaps of map state in the epoch
        pylab.imshow(map[:,:,0])
        pylab.savefig(f"map{e+1}")
        pylab.clf()

        #saving data for animation
        maps[:,:,e] = map[:,:,0]



        if sum(sum(map[:,:,0])) == (N**2)*2: #checking whether everyone is cured to end simulation
            break
    pylab.plot(count[0], label='healthy')
    pylab.plot(count[1], label='infected')
    pylab.plot(count[2], label='cured')
    pylab.legend(loc='upper right')
    pylab.xlabel('epoch')
    pylab.savefig(f"plot.png")
    pylab.clf()

    #preparing for animation
    fig = pylab.figure()
    im = pylab.imshow(maps[:, :, 0])

    def init():
        im.set_data(np.zeros((N, N)))

    def animate(i):
        data = maps[:, :, i]
        im.set_data(data)
        return im

    #animation
    anim = FuncAnimation(fig, animate, init_func=init, frames=M, repeat=False)
    anim.save('spreading.gif', writer='imagemagick')


"""
Była próba wykorzystania biblioteki argparse jednak z poziomu terminala wykrywało dziwne błędy w kodzie, których normalnie nie było + nie widziało biblioteki numpy?
Możliwe, że wyhashowany kod działa, ale nie na moim komputerze, więc wykorzystałam niepreferowane rozwiązanie
"""

#parser = argparse.ArgumentParser()
#parser.add_argument("N", help="size of the map",type=int)
#parser.add_argument("k", help="number of infected cells",type=int)
#parser.add_argument("p_w", help="probability of curing the infection",type=float)
#parser.add_argument("p_z", help="probability of spreading the infection",type=float)
#parser.add_argument("M", help="number of epochs",type=int)
#parser.add_argument("r", help="radius of spreadage",type=int)
#args = parser.parse_args()

#main(args)

#Getting the data for simulation from the user
N = int(input("Set the size of the map (N): "))
k = int(input("Set the number of infected cells (k): "))
p_w = float(input("Set the probability of curing infection (p_w): "))
p_z = float(input("Set the probability of getting infected (p_z): "))
M = int(input("Set how many epochs should the simulation take (M): "))
r = input("Set the radius of spreading the infection (r), if not provided: r=1: ")

if r =='':
    main(N,k,p_w,p_z,M)
else:
    main(N,k,p_w,p_z,M,int(r))
