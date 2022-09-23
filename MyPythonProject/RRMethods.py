import copy

import numpy as np
from numpy import sqrt, dot, cross
from numpy.linalg import norm
import matplotlib.tri as mtri #https://matplotlib.org/stable/gallery/mplot3d/trisurf3d_2.html

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

# https://stackoverflow.com/questions/35176451/python-code-to-calculate-angle-between-three-point-using-their-3d-coordinates
def checkAngle(a, b, c, localAngle):

    # https://thispointer.com/6-ways-to-check-if-all-values-in-numpy-array-are-zero-in-both-1d-2d-arrays-python/
    if (np.all((a == -1))) or (np.all((b == -1))) or (np.all((c == -1))):
        print('Failed, a point being used is one that previously failed; target angle ', localAngle)
        return(0)

    ba = a - b
    bc = c - b

    # print(np.dot(ba, bc))
    # print(np.linalg.norm(ba))
    # print(np.linalg.norm(bc))
    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    if cosine_angle > 1:
        print('Failed, cosine angle >1 (', cosine_angle,') and thus arccos does not exist. Target angle, ',localAngle)
        return(0)

    angle = np.degrees(np.arccos(cosine_angle))

    if(angle > localAngle):
        print('Failed, angle is ', angle, 'which is greater than target angle ', localAngle)
        return(0)
    elif(angle == 0):
        print('Failed, angle is ', angle, ', target angle is', localAngle)
        return(0)
    elif(angle < 1):
        print('Failed, angle is ', angle, 'which although is < target angle ', localAngle,'is too small for the ~ flat surface')
        return(0)
    elif(angle < localAngle):
        print('Pass, angle is ', angle, 'which is less than target angle ', localAngle)
        return(1)
    elif (angle == localAngle):
        print('Pass, angle is equal to ', angle, 'which is equal to target angle ', localAngle)
        return (1)


# Find the intersection of three spheres
# https://stackoverflow.com/questions/1406375/finding-intersection-points-between-3-spheres
# P1,P2,P3 are the centers, r1,r2,r3 are the radii
# Implementation based on Wikipedia Trilateration article.
def trilaterate(P1, P2, P3, r1, r2, r3, condition, angle):

    if(checkAngle(P1,P2,P3,angle) == 0):
        return np.array([-1, -1, -1])

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
        print('temp4 = ', temp4)
        raise Exception("The three spheres do not intersect!");
    z = sqrt(temp4)

    # 2 points of interest
    if condition == 'RHR':
        p_12 = P1 + x * e_x + y * e_y + z * e_z # RHR
    else:
        p_12 = P1 + x * e_x + y * e_y - z * e_z  # LHR

    return p_12



def solveRR(localFlat, N, ax, xindex, yindex, SR3):
    xindexCopy = xindex
    yindexCopy = yindex # duplicates for point I

    error = False

    for x in range(0,N):

        # points A-E are part of the 5 points in the unit 'zigzag' shape
        pointA = np.array([localFlat[0][xindex,yindex],
                  localFlat[1][xindex,yindex],
                  localFlat[2][xindex,yindex]])

        pointB = np.array([localFlat[0][xindex,yindex+1],
                  localFlat[1][xindex,yindex+1],
                  localFlat[2][xindex,yindex+1]])

        pointC = np.array([localFlat[0][xindex-1,yindex+1],
                  localFlat[1][xindex-1,yindex+1],
                  localFlat[2][xindex-1,yindex+1]])

        pointD = np.array([localFlat[0][xindex-1,yindex+2],
                  localFlat[1][xindex-1,yindex+2],
                  localFlat[2][xindex-1,yindex+2]])

        pointE = np.array([localFlat[0][xindex-2,yindex+2],
                  localFlat[1][xindex-2,yindex+2],
                  localFlat[2][xindex-2,yindex+2]])

        # points F,G are 1A points
        try:
            pointF = trilaterate(pointA, pointB, pointC, 2, SR3, SR3, 'RHR', 150)
        except:
            print('pointF failed (try/except)')
            error = True
        else:
            if pointF[0] > 0:
                localFlat = add2Flat(localFlat, pointF, [xindex - 1, yindex])



        try:
            pointG = trilaterate(pointC, pointD, pointE, SR3, SR3, 2, 'RHR', 150)
        except:
            print('pointG failed (try/except)')
            error = True
        else:
            if pointF[0] > 0: #some weird error where F or G = [nan nan nan] and it plots random lines on my graph. this prevents that
                if pointG[0] > 0:
                    localFlat = add2Flat(localFlat, pointG, [xindex - 2, yindex + 1])
                    plotA(pointA, pointB, pointC, pointD, pointE, pointF, pointG, ax)



        # pointH is a 1B point
        try:
            pointH = trilaterate(pointF, pointC, pointG, 1, 2, 1, 'LHR', 60)
        except:
            print('pointH failed (try/except)')
            error = True
        else:
            if pointH[0] > 0: #some weird error where H = [nan nan nan] and it plots random lines on my graph. this prevents that
                localFlat = add2Flat(localFlat, pointH, [xindex - 2, yindex])
                plotB(pointC, pointF, pointG, pointH, ax)

        xindex = xindex - 2
        yindex = yindex +2

    # point I's
    if N > 1: # N > 1
        xindex = xindexCopy #resetting to original values
        yindex = yindexCopy

        for x in range (0,N-1):
            # refind point G
            pointG = np.array([localFlat[0][xindex - 2, yindex + 1],
                               localFlat[1][xindex - 2, yindex + 1],
                               localFlat[2][xindex - 2, yindex + 1]])

            # refind point A/E end of one unit length, start of another
            pointAE = np.array([localFlat[0][xindex - 2, yindex + 2],
                                localFlat[1][xindex - 2, yindex + 2],
                                localFlat[2][xindex - 2, yindex + 2]])

            # refind point F
            pointF = np.array([localFlat[0][xindex - 3, yindex + 2],
                               localFlat[1][xindex - 3, yindex + 2],
                               localFlat[2][xindex - 3, yindex + 2]])

            try:
                pointI = trilaterate(pointG, pointAE, pointF, SR3, 1, SR3, 'RHR', 120)
            except:
                print('pointI failed (try/except)')
                error = True
            else:
                if pointI[0] > 0: #some weird error where I = [nan nan nan] and it plots random lines on my graph. this prevents that
                    localFlat = add2Flat(localFlat, pointI, [xindex - 3, yindex + 1])
                    plotB(pointAE, pointG, pointF, pointI, ax)

            xindex = xindex - 2
            yindex = yindex + 2




        # # for testing purposes only; printing the flat arrays 1 at a time tor easy viewability. values approximated with np.around
        # localFlatCopy = copy.deepcopy(localFlat)
        # for a in range(0, 3):
        #     localFlatCopy[a] = np.around(localFlatCopy[a], decimals=1)
        #     print(localFlatCopy[a])
        #     print('')

        if (error == True):
            print('Trilateration error in the previous iteration ( N =', N, '); recursion stopped prematurely')
        else:
            print('begin recursion')
            solveRR(localFlat, N - 1, ax, 2 * (N - 1), 0, SR3)


        # solveRR(localFlat, N - 1, ax, 2 * (N - 1), 0, SR3)

    return localFlat




def solveRR2(localFlat, N, ax, xindex, yindex, SR3):
    xindexCopy = xindex
    yindexCopy = yindex # duplicates for point I

    error = False

    for x in range(0,N):
        pointB = np.array([localFlat[0][xindex,yindex+1],
                  localFlat[1][xindex,yindex+1],
                  localFlat[2][xindex,yindex+1]])

        pointC = np.array([localFlat[0][xindex-1,yindex+1],
                  localFlat[1][xindex-1,yindex+1],
                  localFlat[2][xindex-1,yindex+1]])

        pointD = np.array([localFlat[0][xindex-1,yindex+2],
                  localFlat[1][xindex-1,yindex+2],
                  localFlat[2][xindex-1,yindex+2]])

        # points F,G are 1A points
        try:
            pointJ = trilaterate(pointB, pointC, pointD, 2, 1, 2, 'LHR', 180)
        except:
            print('pointJ failed (try/except)')
            error = True
        else:
            # if pointJ[0] > 0:
            localFlat = add2Flat(localFlat, pointJ, [xindex, yindex + 2])
            plotB(pointJ, pointB, pointD, pointC, ax)

        xindex = xindex - 2
        yindex = yindex + 2

    if N > 1:
        xindex = xindexCopy
        yindex = yindexCopy

        for x in range(0, N-1):
            pointD = np.array([localFlat[0][xindex - 1, yindex + 2],
                               localFlat[1][xindex - 1, yindex + 2],
                               localFlat[2][xindex - 1, yindex + 2]])

            pointE = np.array([localFlat[0][xindex-2,yindex+2],
                      localFlat[1][xindex-2,yindex+2],
                      localFlat[2][xindex-2,yindex+2]])

            pointB = np.array([localFlat[0][xindex - 2, yindex + 3],
                               localFlat[1][xindex - 2, yindex + 3],
                               localFlat[2][xindex - 2, yindex + 3]])

            pointJLeft = np.array([localFlat[0][xindex, yindex + 2],
                               localFlat[1][xindex, yindex + 2],
                               localFlat[2][xindex, yindex + 2]])

            pointJRight = np.array([localFlat[0][xindex - 2, yindex + 4],
                               localFlat[1][xindex - 2, yindex + 4],
                               localFlat[2][xindex - 2, yindex + 4]])

            try:
                pointK = trilaterate(pointD, pointE, pointB, SR3, 2, SR3, 'RHR', 120)
            except:
                print('pointK failed (try/except)')
                error = True
            else:
                # if pointK[0] > 0:
                localFlat = add2Flat(localFlat, pointK, [xindex - 1, yindex + 3])
                plotB(pointK, pointD, pointB, pointE, ax)

            try:
                pointL = trilaterate(pointJLeft, pointD, pointK, 1, SR3, SR3, 'LHR', 90)
            except:
                print('pointL failed (try/except)')
                error = True
            else:
                # if pointK[0] > 0:
                localFlat = add2Flat(localFlat, pointL, [xindex, yindex + 3])
                plotB(pointL, pointJLeft, pointK, pointD, ax)

            try:
                pointM = trilaterate(pointK, pointB, pointJRight, SR3, SR3, 1, 'LHR', 90)
            except:
                print('pointM failed (try/except)')
                error = True
            else:
                # if pointK[0] > 0:
                localFlat = add2Flat(localFlat, pointM, [xindex-1, yindex + 4])
                plotB(pointM, pointK, pointJRight, pointB, ax)

            xindex = xindex - 2
            yindex = yindex + 2

        if (error == True):
            print('Trilateration error in the previous iteration ( N =', N, '); recursion stopped prematurely')
        else:
            print('begin RR2 recursion')
            solveRR2(localFlat, N - 1, ax, xindexCopy, yindexCopy + 2, SR3)

    return localFlat















def plotA(lPointA, lPointB, lPointC, lPointD, lPointE, lPointF, lPointG, ax):
    xy = [[lPointA[0], lPointA[1]],
          [lPointB[0], lPointB[1]],
          [lPointC[0], lPointC[1]],
          [lPointD[0], lPointD[1]],
          [lPointE[0], lPointE[1]],
          [lPointF[0], lPointF[1]],
          [lPointG[0], lPointG[1]]]
    xy = np.array(xy)

    # manual defining the orientation (connectivity) of the triangles to be plotted (NPQAB :: 01234)
    triangles = [[0, 1, 5],
                 [1, 2, 5],
                 [2, 3, 6],
                 [3, 4, 6]]

    # https://stackoverflow.com/questions/45243563/creating-a-triangulation-for-use-in-matplotlibs-plot-trisurf-with-matplotlib-tr
    triang = mtri.Triangulation(xy[:, 0], xy[:, 1], triangles=triangles)
    z = [lPointA[2], lPointB[2], lPointC[2], lPointD[2], lPointE[2], lPointF[2], lPointG[2]]
    ax.plot_trisurf(triang, z, color=(1, 0, 1), edgecolor='black')



def plotB(lPointC, lPointF, lPointG, lPointH, ax):
    xy = [[lPointC[0], lPointC[1]],
          [lPointF[0], lPointF[1]],
          [lPointG[0], lPointG[1]],
          [lPointH[0], lPointH[1]]]

    xy = np.array(xy)

    # manual defining the orientation (connectivity) of the triangles to be plotted (NPQAB :: 01234)
    triangles = [[0, 1, 3],
                 [0, 2, 3]]

    # https://stackoverflow.com/questions/45243563/creating-a-triangulation-for-use-in-matplotlibs-plot-trisurf-with-matplotlib-tr
    triang = mtri.Triangulation(xy[:, 0], xy[:, 1], triangles=triangles)
    z = [lPointC[2], lPointF[2], lPointG[2], lPointH[2]]
    ax.plot_trisurf(triang, z, color=(1, 0, 1), edgecolor='black')

