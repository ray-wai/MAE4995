import numpy as np
from numpy import sqrt, dot, cross
from numpy.linalg import norm
import matplotlib.tri as mtri #https://matplotlib.org/stable/gallery/mplot3d/trisurf3d_2.html

# Find the intersection of three spheres
# https://stackoverflow.com/questions/1406375/finding-intersection-points-between-3-spheres
# P1,P2,P3 are the centers, r1,r2,r3 are the radii
# Implementation based on Wikipedia Trilateration article.
def trilaterate(P1, P2, P3, r1, r2, r3):
    temp1 = P2 - P1
    e_x = temp1 / norm(temp1)
    temp2 = P3 - P1
    i = dot(e_x, temp2)
    temp3 = temp2 - i * e_x
    e_y = temp3 / norm(temp3)
    e_z = cross(e_x, e_y) #contains info for x and y; not exclusively z (see pyramidPlot2)
    d = norm(P2 - P1)
    j = dot(e_y, temp2)
    x = (r1 * r1 - r2 * r2 + d * d) / (2 * d)
    y = (r1 * r1 - r3 * r3 - 2 * i * x + i * i + j * j) / (2 * j)
    temp4 = r1 * r1 - x * x - y * y
    if temp4 < 0:
        raise Exception("The three spheres do not intersect!");
    z = sqrt(temp4)

    # 2 points of interest
    p_12_a = P1 + x * e_x + y * e_y + z * e_z
    p_12_b = P1 + x * e_x + y * e_y - z * e_z  # theres probably a way to detemrine which point is correct in trilarterate; more research needed. probably with cross and dot functions

    return p_12_a, p_12_b



def solveZigzag(zz, N, ax, localFlat, mode):
    n = (int((len(zz)-1)/2)) # little n is the size of the diagonal strip of pyramids to be made in a single call; n = 3 means 3 pyramids

    for count in range(0, (2 * n) - 1, 2): #each for loop iteration is a single pyramid

        # a is the tip of the pyramid, solved for using 3 "zigzag" points
        if mode == 'left':
            a = trilaterate(zz[count], zz[count + 1], zz[count + 2], 1, 1, 1)
        else:
            a = trilaterate(zz[count + 2], zz[count + 1], zz[count], 1, 1, 1)  # changed order of zz[count+2] and zz[count] so trilaterate cross fnctn works better/correctly

        # a0 = a[0]
        # a1 = a[1]
        #
        #
        # # test to determine which of the two points (a0,a1) distance 1 away from the 3 zigzag points is the one of interest
        # if (a0[2] < 0):
        #     afinal = a1
        # else:
        #     afinal = a0

        afinal = a[0]

        # b is the missing planar point of the pyramid, solved for using opposite corners and the tip of the pyramid
        b = trilaterate(zz[count], afinal, zz[count + 2], 1, 1, 1)
        b0 = b[0]
        b1 = b[1]

        # test to determine which of the two points (b0,b1) distance 1 away from opposite corners and tip is the one of interest
        # we have to calculate distance instead of comparing equality between points because the trilaterate
        # method sometimes returns points such as (0, 0, 0.0000001*E-16) != (0,0,0)
        xdist = (b0[0] - zz[count + 1][0]) ** 2
        ydist = (b0[1] - zz[count + 1][1]) ** 2  # python power operator (squared)
        zdist = (b0[2] - zz[count + 1][2]) ** 2
        dist = sqrt(xdist + ydist + zdist)

        if (dist < 0.1):  # if true, b0 is a point we already have and b1 is of interest (fairly certain bfinal = b1 for RHS)
            bfinal = b1
        else:  # b1 is a point we already have and b0 is of interest (fairly certain bfinal = b0 for LHS)
            bfinal = b0

        #adding the newly found points to the matrix flat array. afinal is tip, bfinal if flat point
        if mode == 'left':
            localFlat = add2Tip(localFlat, afinal, [int((N - n - 1) + (count / 2)), int(count / 2)])
            localFlat = add2Flat(localFlat, bfinal, [int((N-n)+(count/2)),int(count/2)])
        else:
            localFlat = add2Tip(localFlat, afinal, [int(count / 2), int((N - n - 1) + (count / 2))])
            localFlat = add2Flat(localFlat, bfinal, [int(count / 2), int((N - n) + (count / 2))])

        #plotting a single pyramid
        plot(zz[count], zz[count + 1], zz[count + 2], afinal, bfinal, ax)

        # the following 10 lines construct the zigzag for the next recursion, which is called ~15 lines down
        # first iteration opf for loop
        if (count == 0):
            zztemp = np.array([bfinal])

        # not first iteration of for loop
        else:
            zztemp = np.append(zztemp, [bfinal], axis = 0)

        # NOT last iteration
        if (count != (2 * n) - 2):
            zztemp = np.append(zztemp, [zz[count + 2]], axis = 0)

    # if additional recursion is necessary (square pyramid not fully complete)
    if(n-1>0): #n-1
        if mode == 'left':
            localFlat = solveZigzag(zztemp, N, ax, localFlat, 'left')
        else:
            localFlat = solveZigzag(zztemp, N, ax, localFlat, 'right')

    return localFlat



#add2Flat:
    # adds a 3D point (localPoint) to the point array (localFlat). each points' x,y,z coordinates are all inserted into their
    # respective array at the same index.
def add2Flat(localFlat, localPoint, index):
    # localFlat[0] #x-array
    # localFlat[1] #y-array
    # localFlat[2] #z-array
    #
    # index[0] #x-index for array insertion
    # index[1] #y-index for array insertion
    #
    # localPoint[0] #x-coord for point
    # localPoint[1]  #y-coord for point
    # localPoint[2]  #z-coord for point

    for a in range(0,3): #3 iterations: x, y, and z
        localFlat[a][index[0],index[1]] = localPoint[a]

    return localFlat #returning the modified point array back to the main function

def add2Tip(localFlat, localPoint, index):
    for a in range(0,3):
        localFlat[a+3][index[0],index[1]] = localPoint[a] #a+3 instead of a (add2Flat) so that we write to the tip arrays, not the planar ones

    return localFlat

# plot:
    # the following lines of code plot a singular square pyramid based on the
    # two iterations of trilaterate that were just executed
    # this is necessary as ax.plot_trisurf does not plot the pyramids correctly when all pyramids given at once
    # https://stackoverflow.com/questions/12423601/simplest-way-to-plot-3d-surface-given-3d-points
    # placing x and y points of 3D points into an array (simplifies code ~13 lines down in mtri.Triangulation)
def plot(pointN, pointP, pointQ, afinal, bfinal, ax):
    xy = [[pointN[0], pointN[1]],
          [pointP[0], pointP[1]],
          [pointQ[0], pointQ[1]],
          [afinal[0], afinal[1]],
          [bfinal[0], bfinal[1]]]
    xy = np.array(xy)

    # manual defining the orientation (connectivity) of the triangles to be plotted (NPQAB :: 01234)
    triangles = [[0, 1, 3],
                 [1, 3, 2],
                 [2, 3, 4],
                 [4, 3, 0]]

    # https://stackoverflow.com/questions/45243563/creating-a-triangulation-for-use-in-matplotlibs-plot-trisurf-with-matplotlib-tr
    triang = mtri.Triangulation(xy[:, 0], xy[:, 1], triangles=triangles)
    z = [pointN[2], pointP[2], pointQ[2], afinal[2], bfinal[2]]
    ax.plot_trisurf(triang, z, color=(1, 0, 1), edgecolor='black')