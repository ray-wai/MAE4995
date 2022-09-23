import matplotlib.pyplot as plt
import numpy as np
import RRMethods

# Pattern Size (user input):
# so far, N = 4, alpha = 75, and SR3 = np.sqrt(3) is the best results I have found thus far
N = 4 # number of repeating units
plot = 'Yes'
SR3 = np.sqrt(3) #np.sqrt(3) = 1.73205080757

# figure settings
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')
plt.xlim([0, 4*N])
plt.ylim([-2*N, 2*N])
ax.set_zlim(-2*N, 2*N)  # https://stackoverflow.com/questions/37521910/set-zlim-in-matplotlib-scatter3d
ax.set_xlabel("X")
ax.set_ylabel("Y")  # https://likegeeks.com/3d-plotting-in-python/
ax.set_zlabel("Z")

flatX = np.zeros(((2*N)+1)*((2*N)+1)).reshape((2*N)+1,(2*N)+1)
flatY = np.zeros(((2*N)+1)*((2*N)+1)).reshape((2*N)+1,(2*N)+1)
flatZ = np.zeros(((2*N)+1)*((2*N)+1)).reshape((2*N)+1,(2*N)+1)
cost = 0 # this variable could be used to create a single value cost in future research
flat = [flatX, flatY, flatZ, cost]

# projected angles used in the initial determination of the initial zigzag
# https://math.stackexchange.com/questions/2207665/projecting-an-angle-from-one-plane-to-another-plane
# tan gamma = cos alpha * tan beta
alpha = 45 * np.pi/180 #range from 0 to 90
beta = np.pi/6 # 30 degrees
gamma = np.arctan(np.cos(alpha)*np.tan(beta))

# indexes start in bottom left corner, zigzag towards top right
xindex = 2*N
yindex = 0

# initial zigzag
for x in range (0,N):
    pointB = [SR3/2*np.cos(alpha) + flat[0][xindex,yindex], # x
              -0.5*np.cos(alpha) + flat[1][xindex,yindex], # y
              np.sin(alpha) + flat[2][xindex,yindex]] # z
    flat = RRMethods.add2Flat(flat, pointB, [xindex,yindex+1])

    pointC = [SR3*np.cos(np.pi/6-gamma)+flat[0][xindex,yindex+1],
              -SR3*np.sin(np.pi/6-gamma)+flat[1][xindex,yindex+1],
              flat[2][xindex,yindex+1]]
    flat = RRMethods.add2Flat(flat, pointC, [xindex-1, yindex+1])

    pointD = [SR3*np.cos(np.pi/6-gamma)+flat[0][xindex-1,yindex+1],
              SR3 * np.sin(np.pi / 6 - gamma) + flat[1][xindex-1, yindex + 1],
              flat[2][xindex-1,yindex+1]]
    flat = RRMethods.add2Flat(flat, pointD, [xindex-1, yindex+2])

    pointE = [SR3/2*np.cos(alpha) + flat[0][xindex-1,yindex+2],
              0.5*np.cos(alpha) + flat[1][xindex-1,yindex+2],
              -np.sin(alpha) + flat[2][xindex-1,yindex+2]]
    flat = RRMethods.add2Flat(flat, pointE, [xindex-2, yindex+2])

    xindex = xindex-2
    yindex = yindex+2

# The following ~55 lines of code are an alternate way of creating the initial zigzag using "polar" coordinates
# The implementation allows the user to individually define the XY plane angle (theta) and the angle measured from
# the Z axis (phi) for each of the 4 remaining points in the 5 point unit cell since the original point is at the origin
# the implementation still remains limited in that it is essentially impossible to have the points curve along the X-axis,
# as the for loop ensures that each point E is in a line with all other point E's, for example. the same holds for
# each other point in each unit, A-line, B-line, C-line, D-line, and E-line. Further research to figure out
# how to curve the unit cells along the x-axis would certainly be of interest
# # $$$$$$$$$$$$$$$$
#
# thetaAB = 330 * np.pi/180
# phiAB = 90 * np.pi/180
#
# thetaBC = 0 * np.pi/180
# phiBC = 90 * np.pi/180
#
# thetaCD = 0 * np.pi/180
# phiCD = 90 * np.pi/180
#
# thetaDE = 30 * np.pi/180
# phiDE = 90 * np.pi/180
#
# # indexes start in bottom left corner, zigzag towards top right
# xindex = 2*N
# yindex = 0
#
# # initial zigzag
# for x in range (0,N):
#     pointB = [1*np.sin(phiAB)*np.cos(thetaAB) + flat[0][xindex,yindex], # x
#               1*np.sin(phiAB)*np.sin(thetaAB) + flat[1][xindex,yindex], # y
#               1*np.cos(phiAB) + flat[2][xindex,yindex]] # z
#     flat = RRMethods.add2Flat(flat, pointB, [xindex,yindex+1])
#
#     pointC = [SR3*np.sin(phiBC)*np.cos(thetaBC) + flat[0][xindex,yindex+1],
#               SR3*np.sin(phiBC)*np.sin(thetaBC) + flat[1][xindex,yindex+1],
#               SR3*np.cos(phiBC) + flat[2][xindex,yindex+1]]
#     flat = RRMethods.add2Flat(flat, pointC, [xindex-1, yindex+1])
#
#     pointD = [SR3*np.sin(phiCD)*np.cos(thetaCD) + flat[0][xindex-1,yindex+1],
#               SR3*np.sin(phiCD)*np.sin(thetaCD) + flat[1][xindex-1, yindex + 1],
#               SR3*np.cos(phiCD) + flat[2][xindex-1,yindex+1]]
#     flat = RRMethods.add2Flat(flat, pointD, [xindex-1, yindex+2])
#
#     pointE = [1*np.sin(phiDE)*np.cos(thetaDE) + flat[0][xindex-1,yindex+2],
#               1*np.sin(phiDE)*np.sin(thetaDE) + flat[1][xindex-1,yindex+2],
#               1*np.cos(phiDE) + flat[2][xindex-1,yindex+2]]
#     flat = RRMethods.add2Flat(flat, pointE, [xindex-2, yindex+2])
#
#     xindex = xindex-2
#     yindex = yindex+2
#
# # $$$$$$$$$$$$$$$$$$$$$$$






flat = RRMethods.solveRR(flat, N, ax, 2*N, 0, SR3) # 2N is xindex, 0 is yindex
print('starting rr2')
flat = RRMethods.solveRR2(flat, N, ax, 2*N, 0, SR3)

# for testing purposes only; printing the flat arrays 1 at a time tor easy viewability. values approximated with np.around
for a in range(0,3):
    flat[a] = np.around(flat[a], decimals = 1)
    print(flat[a])
    print('')

# test line that can be used to plot the initial zigzag as points
# ax.scatter(flat[0], flat[1], flat[2], color=(1, 0, 1), edgecolor='black')

if plot == 'Yes':
    plt.show()


# distance point testing
# print('')
# cstar = np.arcsin(((np.sin(alpha))/2))
# angleb = (np.pi) - cstar - alpha
# sideb = (2*np.sin(angleb))/(np.sin(alpha))
# # print(sideb)
#
# pointcstar = [SR3/2*sideb, # x
#               -sideb/2, # y
#               0] # z
#
# print(pointcstar)
# print(flat[0][(2*N)-1,1])
# print(flat[1][(2*N)-1,1])
# print(flat[2][(2*N)-1,1])
#
#

# p1 = np.array([pointcstar[0], pointcstar[1], pointcstar[2]])
# p2 = np.array([flat[0][(2*N)-1,1], flat[1][(2*N)-1,1], flat[2][(2*N)-1,1]])
# print(p1)
# print(p2)
# squared_dist = np.sum((p1-p2)**2, axis=0)
# dist = np.sqrt(squared_dist)
# print('distance = ', dist)



