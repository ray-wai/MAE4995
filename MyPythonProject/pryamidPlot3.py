import matplotlib.pyplot as plt
import numpy as np
import methods

# Pattern Size (user input):
N = 5
plot = 'Yes'

# figure settings
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')
plt.xlim([0, N])
plt.ylim([0, N])
ax.set_zlim(0, N)  # https://stackoverflow.com/questions/37521910/set-zlim-in-matplotlib-scatter3d
ax.set_xlabel("X")
ax.set_ylabel("Y")  # https://likegeeks.com/3d-plotting-in-python/
ax.set_zlabel("Z")

# matrices to contain all the organized points of the mesh
flatX = np.zeros((N+1)*(N+1)).reshape(N+1,N+1)
flatY = np.zeros((N+1)*(N+1)).reshape(N+1,N+1)
flatZ = np.zeros((N+1)*(N+1)).reshape(N+1,N+1)
tipX = np.zeros(N*N).reshape(N,N)
tipY = np.zeros(N*N).reshape(N,N)
tipZ = np.zeros(N*N).reshape(N,N)
flat = [flatX, flatY, flatZ, tipX, tipY, tipZ]

# initial zig zag
zigzag = np.array([[0, 0, 0]])
u = np.array([1,0,0])
theta = np.pi/2*1.1 #curvy pyramids
# theta = np.pi/2 #90 degree normal pyramids
v = np.array([np.cos(theta),np.sin(theta),0])

# https://numpy.org/doc/stable/reference/generated/numpy.append.html
# this for loop creates the zigzag points based on the user defined value of N
for x in range(0, N):
    zigzagPoint0 = [zigzag[-1,:]+u]
    zigzag = np.append(zigzag, zigzagPoint0, axis=0)  # [1,0,0], [2,1,0], [3,2,0]...
    flat = methods.add2Flat(flat, zigzagPoint0[0], [x,x+1])

    zigzagPoint1 = [zigzag[-1,:]+v]
    zigzag = np.append(zigzag, zigzagPoint1, axis=0)  # [1,1,0], [2,2,0], [3,3,0]...
    flat = methods.add2Flat(flat, zigzagPoint1[0], [x+1,x+1])

# ax.scatter(flat[0],flat[1],flat[2], color = (1,0,1), marker='o', linewidths = 5)
# ax.scatter(flat[3],flat[4],flat[5], color = (1,0,1), marker='o', linewidths = 5)
# begin recursive calls, solving the geometry. 'flat' is an array containing all the points in the mesh in matrix form
# modified zigzag does not include first or last point of the zigzag for the second call (RHS)
flat = methods.solveZigzag(zigzag, N+1, ax, flat, 'left')
flat = methods.solveZigzag(zigzag[1:len(zigzag)-1], N+1, ax, flat, 'right')

# for testing purposes only; printing the flat arrays 1 at a time tor easy viewability. values approximated with np.around
for a in range(0,6):
    flat[a] = np.around(flat[a], decimals = 1)
    print(flat[a])
    print('')

if plot == 'Yes':
    plt.show()
