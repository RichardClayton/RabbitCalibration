#!/bin/bash
# S1S2 pacing of cell model using bench


#{{{ process bash command line
if [ $# -ne 5 ]
then
    echo "[ERROR]: call as <script_name> D t_in t_out t_open t_close"
    exit
fi
#}}}


# bench command
bench \
  --imp mMS \
  --imp-par="V_gate=0.1,a_crit=0.1,tau_in=$2,tau_out=$3,tau_open=$4,tau_close=$5" \
  --stim-volt 1.0 \
  --stim-dur 1 \
  --restitute S1S2 \
  --res-file data/restitution_protocol.txt \
  --res-trace data/restout_trace.txt \
  --fout=data/restout


python -c \
"
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('data/restout.txt'); t = data[:,0]; V = data[:,1];
plt.plot(t, V);
plt.axhline(0, linestyle = '--', alpha = 0.5);
plt.show()

"

