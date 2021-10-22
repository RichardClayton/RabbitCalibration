# hdf5Utils

Utilities for handling hdf5 files with python.


## TODO

* display attributes of chosen object after selection?


## Example

It is useful to define a script that can be called from the command line:

```python
#!/usr/bin/env python

import sys
import os

import hdf5utils as hu

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Error: call as ./{} <filename>".format(os.path.basename(__file__)))
    else:
        hu.show_menu(sys.argv[1])
```

