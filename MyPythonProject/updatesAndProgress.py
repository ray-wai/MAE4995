# 11/18/21
    # plotting is good, fixed the error with trilaterate where the mesh was at times inverting.
    # the trilaterate fix is not a rigid one, but it does work with the current pyramid structure
    # main issue is that trilaterate returns 2 points, but only 1 is of interest. right now,
    # just the cross product determines the orientation. This is probably a problem that should be addressed when the
    # points are transitioned into vectors

    # necessary improvements: need to support point output

# 2/11/21
    # point output now supported with newLeftPoionts, newRightPoints, and allPoints. solveZigzag left and right
    # were modified so that the recursive calls correctly append the new points to the 1D array

    # next steps (met with Nassar @10am 2/11/21
    # need to place the points into 6, 2D arrays:
        # 3 of them are for x, y, z, planar points
        # 3 of them are for x, y, z, tips of pyramids
    # discussion necessitated that point organization is necessary before switching to vectors

    # necessary improvements: modify/organize point output s/t it follows the 6, 2D array structure and
    # thus connectivity will be implicitly defined
    # possibly utilize numpy.diag or diagonal

# 2/14/22
    # was able to get the flat array created. zigzag insertion works correctly with new method add2Flat
    # need to continue and get the leftside points added in the recursion. I think modification of the function parameters
    # to include a never changing constant N+1 will help, particularly in that i'm not ussing the diag function.
    # I also need to send the flat array to the method, not the newPoints lists.
    # first point in the left array should have index [(Nstar-n)+count/2, count/2]

# 2/16/22
    # methods line 144, indexing starts at 1, not 0
    # this is a problem rooted in methods line 194
    # what might be most useful is to modify the original zigzag and remove the first (0th) and last point, so that
    # a filler point [0,0,0] is not needed to be used. might be possible to combine solveZigZag left and right into a single
    # that gets called twice; once for the left side and once for the right side, but the initial zigzags would be different (RHS smaller)

# 2/17/22
    # flat array fully implemented with 6, 2D arrays. working correctly.
    # mode was added to solveZigzag so that i could delete solveZigzag right.
        # the modified zigzag input to solveZigzagRight resulted in a lot of the same code being used as solveZigzag(Left) which is why this was done
    # plotting functions in solveZigzag were separated into their own method for cleanliness.
    # add2Tip and add2Flat methods added to add new points to the matrix array

    # things to discuss at next meeting:
        # share progress on matrix output for points
        # need to at least mention the email i got from susie on presenting at the research forum
        # next steps?:
            # energy calculation?
            # vectorize?
            # other meshes? (different shapes other than 1x1 pyramids)

# 2/21/22
    # decided to move on to the ron resh mesh.

# 3/2/22
    # issue with N>1 section not working properly
    # every third set of triangles in the second iteration is angled weird.
        # recursion doesnt really work for this reason
    # not sure if i need to define the initial zigzag differently
    # https://www.youtube.com/watch?v=UXENKmAUL0E&list=RDUXENKmAUL0E&start_radio=1&rv=UXENKmAUL0E&t=101