"""
   h5pyUtils.py

   Module with utility functions for handling HDF5 files via h5py module.

   Created: 10-Feb-2020
   Author:  Sam Coveney

"""


from __future__ import print_function
import os.path
import h5py
import numpy as np


def createHDF5(filename):
    """Create an empty HDF5 file with some groups but no datasets."""

    if os.path.isfile(filename):
        print ("File {:s} exists, so doing nothing.".format(filename))
        return
    else:
        print ("File does not exist, creating new file.".format(filename))

    FILE = h5py.File(filename, "a")

    FILE.close()


def scanHDF5(filename, recursive=True, tab_step=2):
    """Display the contents of a HDF5 file.
    
       Modified version of this: https://stackoverflow.com/questions/43371438/how-to-inspect-h5-file-in-python 
    """
    def scan_node(g, tabs=0):

        if tabs > 0:
            print(' ' * tabs + "|--", g.name.split("/")[-1])
        else:
            print("\n=" + "="*len(filename) + "=" + "\n {:s} ".format(filename) + "\n-" + "-"*len(filename) + "-" )

        for k, v in g.items():
            if isinstance(v, h5py.Dataset):
                print(' ' * tabs + ' ' * tab_step + '  >', v.name.split("/")[-1])
            elif isinstance(v, h5py.Group) and recursive:
                scan_node(v, tabs = tabs + tab_step)

    with h5py.File(filename, 'r') as f:
        scan_node(f)
        f.close()
        print("=" + "="*len(filename) + "=\n")


def deleteGroup(filename, group):
    """Delete a group."""

    FILE = h5py.File(filename, "r+")

    try:
        del FILE[group]
        print("[DELETE]: <{:s}> group deleted.".format(group))
    except:
        pass

    FILE.close()


def createGroup(filename, group, subgroups = []):
    """Create a new group after first deleting the old one."""

    #deleteGroup(filename, group)
    try:
        FILE = h5py.File(filename, "r+")
        GROUP = FILE.create_group(group)
        print("[CREATE]: <{:s}> group created.".format(group))

        for sg in subgroups:
            GROUP.create_group(sg)
    except:
        print("Group exists, so not creating a new one.")
        pass

    FILE.close()


def getDataset(filename, group, dataset):
    """Return a NumPy dataset under the supplied group name."""

    FILE = h5py.File(filename, "r")

    GROUP = FILE[group]

    try:
        data = np.array(GROUP[dataset])
    except:
        print("[ERROR]: <{:s}> dataset in <{:s}> group does not exist, returning None.".format(dataset, group))
        return None

    FILE.close()

    return data


def deleteDataset(filename, group, dataset):
    """Delete a dataset under the supplied group name."""

    FILE = h5py.File(filename, "r+")

    GROUP = FILE[group]

    try:
        del GROUP[dataset]
        print("[DELETE]: <{:s}> dataset in <{:s}> group deleted.".format(dataset, group))
    except:
        pass

    FILE.close()


def createDataset(filename, group, dataset, data):
    """Create a NumPy dataset under the supplied group name."""

    deleteDataset(filename, group, dataset)

    FILE = h5py.File(filename, "r+")

    GROUP = FILE[group]

    GROUP.create_dataset(dataset, data = data)

    print("[CREATE]: <{:s}> dataset in <{:s}> group created.".format(dataset, group))

    FILE.close()


