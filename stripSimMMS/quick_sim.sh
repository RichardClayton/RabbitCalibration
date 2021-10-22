#!/bin/bash

#{{{ process bash command line
if [ $# -ne 8 ]
then
    echo "[ERROR]: call as <script_name> GEOMETRY D t_in t_out t_open t_close S1 num_beats"
    exit
fi
#}}}

if [[ "$1" == "square" ]]; then

    echo "Running square simulation with ${8} x S1: ${7} beats"
    ./carp_square.sh $2 $3 $4 $5 $6 $7 $8; python -c "import run_sim; CV, APD = run_sim.calc_CV_APD('square', $7, $8); print('{:1.3f}  {:3.1f}'.format(CV, APD))"

else

    echo "Running strip simulation with ${8} x S1: ${7} beats"
    ./carp_strip.sh $2 $3 $4 $5 $6 $7 $8; python -c "import run_sim; CV, APD = run_sim.calc_CV_APD('strip', $7, $8); print('{:1.3f}  {:3.1f}'.format(CV, APD))"

fi

