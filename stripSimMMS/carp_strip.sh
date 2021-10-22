#!/bin/bash 
#
# Tissue *strip* simulation to create parameter:APD:CV table
# Author: Sam Coveney
# Date:   27-06-2020
#
# Call this script as "./carp_strip.sh D t_in t_out t_open t_close S1"
#
# The simulation mesh can be created with "./make_strip.sh"
# NOTE: the strip coordinates are 'hard-coded' into *this* script by means of stimulus specification.
#


#echo -e "============================"
#echo "CARP tissue strip simulation"
#echo -e "============================\n"


#{{{ process bash command line
if [ $# -lt 7 ]
then
    echo "[ERROR]: call as <script_name> D t_in t_out t_open t_close S1 num_beats (S2)"
    exit
fi

# read in S1 value from command line
S1=${6}
NUM_STIM=${7}

if [ $# -ge 8 ]
then
    DYN=1
    S2FIRST=$9
    S2LAST=${10}
    #echo "Using dynamic S2 protocol"
else
    DYN=0
    S2FIRST=-1
    S2LAST=-1
fi
#}}}


#{{{ Simulation parameters
#############################

# seems to give same answers (when prepacing) as 50, 10
TIMESTEP=100  # time step (microseconds)
ODESUB=5    # nb of sub-iterations for ionic current;

## OUTPUT ##
#SPACEDT=${TEND}   # writes the solution in a file every x milliseconds
#SPACEDT=5   # writes the solution in a file every x milliseconds
TIMEDT=50    # writes screen output every x milliseconds
TSHACT=0.7  # Threshold to evaluate DEPOLARISATION 
TSHREP=0.2  # Threshold to evaluate REPOLARISATION

## STIMULUS ##
ISTIM=2.0       # intensity
TDUR=1.0        # duration in ms

# all the same stimuli parameters
NPLSS1=1        # nb of stimuli applied to the 1st stimulus
BCLS1=1       # period for stimulus 1 # NOTE: since I only use 1 beat per defined stimulus, I set this to 1

#}}}


#{{{ mesh, param files, and output dirs
#######################################

PBS_O_WORKDIR=$PWD

MESHNAME='strip'
MESHDIR=${PBS_O_WORKDIR}'/'${MESHNAME}
MESH=${MESHDIR}'/'${MESHNAME}
PARDIR=${PBS_O_WORKDIR}'/'${MESHNAME}
#PAR=${PARDIR}'/square.par'
#}}}


#{{{ settings Desktop vs HPC
#############################

# run on Desktop 
#echo "Desktop settings"
PBS_O_WORKDIR=$PWD
NPROC=1
CARP='carpentry'
LAUNCHER='mpirun'

CARPCOMMOSTRING="${LAUNCHER} -n ${NPROC} ${CARP} -meshname ${MESH} -dt ${TIMESTEP} -ode_fac ${ODESUB} -mass_lumping 0 -bidomain 0 -bidm_eqv_mono 1 -parab_solve 1 -cg_precond 2 -cg_tol_parab 1.0e-7"

# change to output directory
cd $PBS_O_WORKDIR
#}}}

    
#{{{ Conductivity and Cell parameters
#####################################
#echo "Setting parameters"
#echo "------------------"

DIFFUSION=${1}
TAU_IN=${2}
TAU_OUT=${3}
TAU_OPEN=${4}
TAU_CLOSE=${5}
 
#echo "Diffusion: ${DIFFUSION}"
#echo "Tau_in: ${TAU_IN}"
#echo "Tau_out: ${TAU_OUT}"
#echo "Tau_open: ${TAU_OPEN}"
#echo "Tau_closed: ${TAU_CLOSE}"
#echo ""

#}}}


#{{{ DEFINE THE STIMULUS

TEND=$(((${NUM_STIM}+1)*$S1))     # end time (milliseconds) # FIXME: num_stim * S1 is probably okay for strip, although not for atrium
SPACEDT=${TEND}   # writes the solution in a file every x milliseconds

ONAME=${PBS_O_WORKDIR}'/data/output_'${S1}

#{{{ define location of stimulus
STRIP_START_X=-300.0
STRIP_WIDTH=600.0
#STRIP_START_Y=-6000.0
#STRIP_START_Y=-24000.0
# automatically determine from .pts file
read -ra ARR  <<< `head -2 strip/strip.pts | tail -1`
STRIP_START_Y=${ARR[1]}
#echo "Strip begins at: $STRIP_START_Y"
#echo ${STRIP_START_Y}
#}}}

TOTAL_STIMULUS=""

if [ $S1 -ge 350 ]; then
#{{{ original stimulus of all S1 beats

    #echo "Using standard S1 protocol"

    for i in $(seq 0 $((${NUM_STIM}-1))); do

      S1START=$(($S1*$i))

      STIM_LOC="-stimulus[$i].x0 $STRIP_START_X -stimulus[$i].xd $STRIP_WIDTH -stimulus[$i].y0 $STRIP_START_Y -stimulus[$i].yd 500."

      STIM_COMMON="-stimulus[$i].stimtype 0 -stimulus[$i].strength ${ISTIM} -stimulus[$i].duration ${TDUR} -stimulus[$i].npls ${NPLSS1} -stimulus[$i].bcl ${BCLS1} -stimulus[$i].dump_vtx_file 0 -stimulus[$i].start ${S1START} $STIM_LOC"

      TOTAL_STIMULUS=${TOTAL_STIMULUS}" $STIM_COMMON $STIM_LOC -stimulus[$i].name S${i} " # NOTE: needs whitespace on end

      S1_PLUS=0 # for this standard S1 protocol, set to zero
    done
#}}}
else
    #{{{ new stimulus idea, starting with larger S1 and reducing it down

    #echo "Using dynamic S1 reduction protocol"

    S1_PLUS=50 # how much higher above S1 do we begin?
    S1_PRE=2 # how many larger 'easing in' S1s do we want?
    FIRST_BEATS=1 # how many S1 + S1_PLUS beats to do before decrease begins?

    for i in $(seq 0 $((${NUM_STIM}+${S1_PRE}+${FIRST_BEATS}))); do

      if [ $i -le $(( 0 + ${FIRST_BEATS} )) ]; then
          S1START=$(( ($S1 + $S1_PLUS)*$i ))
          #if [ $i -eq 0 ]; then 
          #  echo "S1START (first beat): $S1START"
          #else
          #  echo "S1START (pre beats): $S1START"
          #fi
      else

          TMP=$(( ($S1 + $S1_PLUS  - ($S1_PLUS/$S1_PRE)*($i-1 - ${FIRST_BEATS})  ) ))
          if [ $TMP -le $S1 ]; then
            TMP=$S1;
            S1START=$(( $S1START + $TMP ))
            #echo "S1START (final S1): $S1START"
          else
            S1START=$(( $S1START + $TMP ))
            #echo "S1START (decreasing S1): $S1START"
          fi
      fi

      STIM_LOC="-stimulus[$i].x0 $STRIP_START_X -stimulus[$i].xd $STRIP_WIDTH -stimulus[$i].y0 $STRIP_START_Y -stimulus[$i].yd 500."

      STIM_COMMON="-stimulus[$i].stimtype 0 -stimulus[$i].strength ${ISTIM} -stimulus[$i].duration ${TDUR} -stimulus[$i].npls ${NPLSS1} -stimulus[$i].bcl ${BCLS1} -stimulus[$i].dump_vtx_file 0 -stimulus[$i].start ${S1START} $STIM_LOC"

      TOTAL_STIMULUS=${TOTAL_STIMULUS}" $STIM_COMMON $STIM_LOC -stimulus[$i].name S${i} " # NOTE: needs whitespace on end

    done

    # update TEND
    #TEND=$(($S1START + $S1))
    TEND=$(($S1START + 400)) # NOTE: use 1000 when strip is extremely long
    SPACEDT=$TEND
    #echo "TEND: $TEND"

    # must update number of stimuli 
    NUM_STIM=$((${NUM_STIM}+${S1_PRE}+${FIRST_BEATS}+1)) 
    #echo "NUM_STIM: ${NUM_STIM}"

    #}}}
fi

#}}}


#{{{ define the S2 stimulus

    TSAVE=$(($S1START + 100))

    if [ $DYN -eq 1 ]
    then
        STIM_LOC="-stimulus[0].x0 $STRIP_START_X -stimulus[0].xd $STRIP_WIDTH -stimulus[0].y0 $STRIP_START_Y -stimulus[0].yd 500."

        STIM_COMMON="-stimulus[0].stimtype 0 -stimulus[0].strength ${ISTIM} -stimulus[0].duration ${TDUR} -stimulus[0].npls ${NPLSS1} -stimulus[0].bcl ${BCLS1} -stimulus[0].dump_vtx_file 0 $STIM_LOC"

    fi

#}}}


#{{{ create parameter command line strings

IMP_GENERAL="-num_imp_regions 1 -imp_region[0].name MYO -imp_region[0].cellSurfVolRatio 1.0 -imp_region[0].im mMS"
IMP_PARAM="-imp_region[0].im_param V_gate=0.1,a_crit=0.1,tau_in=$TAU_IN,tau_out=$TAU_OUT,tau_open=$TAU_OPEN,tau_close=$TAU_CLOSE"
G_GENERAL="-num_gregions 1 -gregion[0].name MYO_0 -gregion[0].num_IDs 1 -gregion[0].ID[0] 1"

# when bidomain equivalent, we need to multiply the diffusion by 2.0
G_VAL=`python -c "print('{0:0.4f}'.format(2.0*$DIFFUSION*1000.0))"`
#echo "G_VAL is $G_VAL"

G_PARAM="-gregion[0].g_el $G_VAL -gregion[0].g_et $G_VAL -gregion[0].g_il $G_VAL -gregion[0].g_it $G_VAL"
PARAMS="$IMP_GENERAL $IMP_PARAM $G_GENERAL $G_PARAM"

#}}}


#{{{ define LAT recordings
LATSSTRING1="-lats[0].measurand 0 -lats[0].all 1 -lats[0].method 1 -lats[0].mode 0 -lats[0].threshold ${TSHACT} -lats[0].ID tact_${TSHACT}"
LATSSTRING2="-lats[1].measurand 0 -lats[1].all 1 -lats[1].method 1 -lats[1].mode 1 -lats[1].threshold ${TSHREP} -lats[1].ID trep_${TSHREP}"
#}}}


#{{{ prepacing cell model, same states for all tissue; NOTE: does not run tissue simulation, but appends initialize to the CARP command
if false; then
    BENCH_S1=$(( $S1 + 50 ))
    STIM_DUR_BENCH=1.0
    CELL_NUM_STIM=10
    CELL_IMP=" --imp=mMS --imp-par=V_gate=0.1,a_crit=0.1,tau_in=$TAU_IN,tau_out=$TAU_OUT,tau_open=$TAU_OPEN,tau_close=$TAU_CLOSE"
    CELL_STIM="--bcl ${BENCH_S1} --numstim ${CELL_NUM_STIM} --duration $((${CELL_NUM_STIM}*${BENCH_S1})) --stim-volt ${ISTIM} --stim-dur ${STIM_DUR_BENCH}"

    PREPACE="bench $CELL_IMP $CELL_STIM --save-ini-file data/INIT.sv --save-ini-time $((${NUM_STIM}*${BENCH_S1})) --dt 0.10 --dt-out 0.10 --fout"
    echo $PREPACE
    ${PREPACE} # run bench simulation

    # prepacing cell states to be used for tissue
    #TOTAL_STIMULUS="${TOTAL_STIMULUS} -imp_region[0].im_sv_init data/INIT.sv"

#echo "[WARNING]: going to exit"
#exit 1
fi
#}}}


#{{{ run CARP simulation
#if [ $S1 -lt 350 ]; then
if false; then

    #{{{ prepacing cell model, LAT-based states for tissue

    #{{{ single beat simulation to get LAT sequence
    #echo "Single beat simulation to get LAT sequence"

    # FIXME: what if this LAT sequence is also part of the problem? What if we prepaced the state first using the non-LAT-based method?


    # define a single beat, to give us LAT
    S1START=0
    STIM_LOC="-stimulus[0].x0 $STRIP_START_X -stimulus[0].xd $STRIP_WIDTH -stimulus[0].y0 $STRIP_START_Y -stimulus[0].yd 500."
    STIM_COMMON="-stimulus[0].stimtype 0 -stimulus[0].strength ${ISTIM} -stimulus[0].duration ${TDUR} \
                 -stimulus[0].npls ${NPLSS1} -stimulus[0].bcl ${BCLS1} -stimulus[0].dump_vtx_file 0 -stimulus[0].start ${S1START} $STIM_LOC"
    SINGLE_STIMULUS=" $STIM_COMMON $STIM_LOC -stimulus[0].name S0 " # NOTE: needs whitespace on end

    # set different LAT recordings to get single S1 only
    SINGLE_LATSSTRING1="-lats[0].measurand 0 -lats[0].all 0 -lats[0].method 1 -lats[0].mode 0 -lats[0].threshold ${TSHACT} -lats[0].ID tact_${TSHACT}"
    SINGLE_LATSSTRING2="-lats[1].measurand 0 -lats[1].all 0 -lats[1].method 1 -lats[1].mode 1 -lats[1].threshold ${TSHREP} -lats[1].ID trep_${TSHREP}"

    # a single beat, to get the LAT sequence; NOTE: to get the wave to cross the strip in 1 beat, we don't need much time
    ${CARPCOMMOSTRING} -simID ${ONAME} ${PARAMS} -tend 50.0 -spacedt 50.0 -timedt ${TIMEDT} -num_LATs 2 ${SINGLE_LATSSTRING1} ${SINGLE_LATSSTRING2} \
                       -num_stim 1 ${SINGLE_STIMULUS} -dt 10 -ode_fac 1  -imp_region[0].im_sv_init data/INIT.sv  #> /dev/null 2>&1
    #}}}

    #{{{ run prepaced CARP simulation with small dt for short time, save state

    PREPACE_BCL=$(($S1 + ${S1_PLUS})) # NOTE: let's use a large S1 for the prepacing of the cell model
    PREPACE_BEATS=10
    
    SINGLE_STIMULUS="${SINGLE_STIMULUS} -prepacing_beats ${PREPACE_BEATS} -prepacing_lats ${ONAME}/init_acts_tact_0.7-thresh.dat -prepacing_bcl ${PREPACE_BCL}"
    STIM_DUR_PREPACE=1.0
    SINGLE_STIMULUS=" ${SINGLE_STIMULUS} -stimulus[0].duration ${STIM_DUR_PREPACE} -stimulus[0].bcl ${STIM_DUR_PREPACE} "
    PRE_SIM_TIME=1
    PRE_SIM_DT=0.11

    #echo "Prepaced short simulation with small timestep (${PREPACE_BEATS} beats, dt: ${PRE_SIM_DT}, BCL: ${PREPACE_BCL})"
    ${CARPCOMMOSTRING} -simID ${ONAME} ${PARAMS} -tend ${PRE_SIM_TIME} -spacedt ${PRE_SIM_TIME} -timedt ${PRE_SIM_TIME} -num_LATs 2 ${LATSSTRING1} ${LATSSTRING2} \
                       -num_stim 1 ${SINGLE_STIMULUS} -dt ${PRE_SIM_DT} -ode_fac 1 \
                       -num_tsav 1 -tsav[0] ${PRE_SIM_TIME} -write_statef ${ONAME}"/statefile"  #> /dev/null 2>&1
    #}}}

    #{{{ continue CARP simulation with regular dt
    #echo "The rest of the CARP simulation (${NUM_STIM} beats)"

    ${CARPCOMMOSTRING} -simID ${ONAME} ${PARAMS} -tend ${TEND}.0 -spacedt ${SPACEDT} -timedt ${TIMEDT} -num_LATs 2 ${LATSSTRING1} ${LATSSTRING2} \
                       -num_stim ${NUM_STIM} ${TOTAL_STIMULUS} \
                       -start_statef ${ONAME}"/statefile.${PRE_SIM_TIME}.0"  #> /dev/null 2>&1
    #}}}

    #}}}

else
    # tissue pacing only

    #echo "NUM_STIM: $NUM_STIM"
    #echo "SPACEDT: $SPACEDT"
    #echo "TEND: $TEND"
    #SPACEDT=5

    #echo ${ONAME}"/statefile"
    #echo $TSAVE

    # S1 beats
    # --------
    # saves the state just after the last S1 beat at time TSAVE, but continues simulation to time TEND
    ${CARPCOMMOSTRING} -simID ${ONAME} ${PARAMS} -tend ${TEND}.0 -spacedt ${SPACEDT} -timedt ${TIMEDT} -num_LATs 2 ${LATSSTRING1} ${LATSSTRING2} \
                       -num_stim ${NUM_STIM} ${TOTAL_STIMULUS} \
                       -num_tsav 1 -tsav[0] ${TSAVE} -write_statef ${ONAME}"/statefile"  #> /dev/null 2>&1

    # FIXME: WHY DOES IT SLOW DOWN AFTER SAVING THE STATE? I NEED TO SAVE IT AT A DIFFERENT TIME TO WHEN IT ENDS...
    
    # S2 beats
    # --------
    if [ $DYN -eq 1 ] # && [ $S1 -eq 600 ]
    then

        # overly basic S2 selection, not related to S1 at all, or to whether S1 was even successful
        #for S2 in 400 350 300 250 200
        for S2 in `seq -s' ' $S2FIRST -10 $S2LAST`
        do

            S2START=$(($S1START + $S2)) # FIXME: will probably move this to a loop below the S1 simulation, but base off command line for now
            S2_STIMULUS="$STIM_COMMON $STIM_LOC -stimulus[0].name S2PREMATURE -stimulus[0].start ${S2START} "
            #echo "S2_STIMULUS: ${S2_STIMULUS}"

            # update TEND since we have another beat to perform
            TEND_S2=$(($TEND + $S2))

            # check values
            #echo $S1START
            #echo $S2_STIMULUS
            #echo $TSAVE
            #echo $TEND

            # because it is a simulation continuation, it seems to append to the LAT file for S2, which I don't want
            rm -r ${ONAME}_${S2}

            # could now loop over different S2 values, saving to a different output directory, so that all S2 is done here before exiting the simulation?
            ${CARPCOMMOSTRING} -simID ${ONAME}_${S2} ${PARAMS} -tend ${TEND_S2}.0 -spacedt ${SPACEDT} -timedt ${TIMEDT} -num_LATs 2 ${LATSSTRING1} ${LATSSTRING2} \
                               -num_stim 1 ${S2_STIMULUS} \
                               -start_statef ${ONAME}"/statefile.${TSAVE}.0"  > /dev/null 2>&1

        done
    fi

fi
#}}}


