# resultsToHDF5.py
#   Store the strip simulation results in a HDF5 file.

import argparse
import glob
import numpy as np

import hdf5utils as hu

def main(args):

    filename = args.filename

    try:
        params = np.loadtxt(filename)
    except OSError as e:
        print("[ERROR]: " + filename + " not found")
        exit()

    head = filename.split(".")[0]

    hdf5file = head + ".hdf5"
    hu.createHDF5(hdf5file)
    hu.createDataset(hdf5file, "/", "params", params)

    res_files = glob.glob(head + "_*")
    print(res_files)

    for rf in res_files:
        print(rf)
        S1 = rf.split("S1")[-1].split("_")[0].split(".")[0]
        S2 = rf.split("S2")[-1].split("_")[0].split(".")[0]
        print("S1:", S1)
        print("S2:", S2)

        data = np.loadtxt(rf)
        hu.createGroup(hdf5file, S1)
        hu.createDataset(hdf5file, S1, S2, data)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = "Store the strip simulation results in a HDF5 file.")

    parser.add_argument("--filename", type = str, required = True, help = "Parameters filename, including path (relative or absolute).")

    args = parser.parse_args()

    main(args)
