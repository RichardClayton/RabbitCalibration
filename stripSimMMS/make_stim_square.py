"""
    make_stim_square.py

    Create a quarter-circle stimulus in the corner of the square.
"""

import numpy as np
from scipy.spatial.distance import cdist


# read in vertices
X = np.loadtxt("square/square.pts", skiprows = 1)

# vertex in one corner
X_corner = X[0]

# calculate distance from corner to all other vertices
dists = cdist(X_corner[None,:], X).flatten()
#print("dists:", dists)

# stimulation vertices should be with R radius of the corner
R = 3000
condition = (dists <= R)
#print("condition:", condition)

stim_verts = np.argwhere(condition).flatten()
print("Apply stimulus to nodes:", stim_verts)


f = open("corner.vtx", "w+")
f.write("{:d}\n".format(stim_verts.shape[0]))
f.write("intra\n".format(stim_verts.shape[0]))
for s in stim_verts:
    f.write("{:d}\n".format(s))
f.close()

#np.savetxt("strip/strip.pts", X_rot, header = str(X.shape[0]), comments = "")

