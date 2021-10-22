from setuptools import setup

setup(name = 'hdf5utils',
      version = '1.0',
      description = 'Utilities for handling hdf5 files in Python.',
      #url = 'http://github.com/samcoveney/quLATi',
      author = 'Sam Coveney',
      author_email = 'coveney.sam@gmail.com',
      license = 'GPL-3.0+',
      packages = ['hdf5utils'],
      install_requires = [
          'future',
          'urwid',
          'h5py'
      ],
      zip_safe = False)

