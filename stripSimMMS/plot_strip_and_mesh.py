"""Plot the strip and the atrial mesh together."""


import hdf5utils as hu
import pyvista as pv
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation
import trimesh


plotter = pv.Plotter()

# plot atrium
filename = "/home/sam/Documents/DATA/UAC_meshes/M2/M2.hdf5"
X = hu.getDataset(filename, "sim_mesh", "X")
X -= X.mean(axis = 0) 
Tri = hu.getDataset(filename, "sim_mesh", "Tri")
plt_surf = pv.PolyData(X, np.hstack([ np.full(Tri.shape[0], 3)[:,None] , Tri ]))
plotter.add_mesh(plt_surf, show_edges = True, opacity = 1.0, color = "white")

mesh = trimesh.Trimesh(X, Tri, process = False)
edge = mesh.edges_unique_length.mean()
edge_std = mesh.edges_unique_length.std()
print("atrial edge length:", edge, "+/-", edge_std)


# plot strip
X_strip = np.loadtxt("strip/strip.pts", skiprows=1)
X_strip -= X_strip.mean(axis = 0) 
#X_strip[:,2] += 50000

R = Rotation.from_euler('z', 90, degrees = True)
X_strip = R.apply(X_strip)
R = Rotation.from_euler('x', 90, degrees = True)
X_strip = R.apply(X_strip)
X_strip[:,1] -= 27500
X_strip[:,2] -= 15000
X_strip[:,0] -= 12000

# rotation matrix
df = pd.read_csv("strip/strip.elem", skiprows=1, header = None, delim_whitespace=True)
Tri_strip = df.values[:,[1,2,3]]

mesh_strip = trimesh.Trimesh(X_strip, Tri_strip, process = False)
edge_strip = mesh_strip.edges_unique_length.mean()
edge_std_strip = mesh_strip.edges_unique_length.std()
print("strip edge length:", edge_strip, "+/-", edge_std_strip)

plt_strip = pv.PolyData(X_strip, np.hstack([ np.full(Tri_strip.shape[0], 3)[:,None] , Tri_strip ]))
plotter.add_mesh(plt_strip, show_edges = True, opacity = 1.0, scalars = X_strip[:,0], cmap = "jet", clim = [X_strip[:,0].min(),X_strip[:,0].max()])


plotter.show()

