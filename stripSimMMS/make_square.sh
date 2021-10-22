#!/bin/bash
# create a square of tissue of triangular elements
# NOTES:
# automatically zero centres
# the mesh units will be 1000x the sizes specified here

SIZE=4.80

mesher -size[0] $SIZE -size[1] $SIZE -size[2] 0.0 -resolution[0] 300.0 -resolution[1] 300.0 -mesh square -tri2D

mkdir -p square/

mv square.* square/.

echo "Running python script to create stimulus file"
python ./make_stim_square.py

