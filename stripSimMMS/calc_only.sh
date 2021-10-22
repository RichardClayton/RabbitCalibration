#!/bin/bash

#{{{ process bash command line
if [ $# -ne 3 ]
then
    echo "[ERROR]: call as <script_name> GEOMETRY S1 num_beats"
    exit
fi
#}}}

if [[ "$1" == "square" ]]; then

    python -c "import run_sim; CV, APD = run_sim.calc_CV_APD('square', $2, $3); print('{:1.3f}  {:3.1f}'.format(CV, APD))"

else

    python -c "import run_sim; CV, APD = run_sim.calc_CV_APD('strip', $2, $3); print('{:1.3f}  {:3.1f}'.format(CV, APD))"

fi

