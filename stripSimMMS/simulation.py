""" simulation.py

    Use python and carputils to launch CARP simulation.

    This simulation can be:
    * a tissue strip paced from the edge (planar wave)
    * a tissue square paced from the corner (circular wave)

    By supplying --S2 VAL on the command line, this code actually does S1S2S3 pacing.
    The S2 restitutions will then really be S3 rests for given S1 and S2.

    NOTES from 19 Jan 2019 onwards
    * modified APD downstroke detector to 0.1 i.e. APD90
    * modified the strength of stimulation from 2.0 to 1.0
    * modifying fit bounds to allow 'd' to be negative
    * modifying fit bounds to allow 'a' to be negative

    TO CHECK:
    * does calc_CV_APD() really ensure that we have a new propagation? -- YES

"""

#{{{ import modules
import os
import sys

from datetime import date
from datetime import datetime
from carputils import settings
from carputils import tools
from carputils.clean import clean_directory, overwrite_check

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import functools
import pickle
from scipy.optimize import curve_fit
#}}}


#{{{ calculate CV and APD
def calc_CV_APD(job_dir, X, Tri, geometry, num_stim, S1, mode = "S1", S2 = 0):
    """ Calculate CV and APD using results of simulation.

        geometry -- used to determine whether a rotation is necessary.
        num_stim -- used to determine if all the activations occurred.
        S1 -- value of S1
        mode -- should be S1 or S2, simply determines the names of LAT files that get loaded
        S2 -- value of S2, should be set to the S2 value if mode == S2

        NOTE: using mode = "S1" and S2 > 0, this routine will work with the S1S2S3 protocol

        (ERP is determined outside of this routine, based on whether activation occured or not)

        This routine works by discarding all activations before the time of the premature beat.
        The assumption, therefore, is that the propagation for the beat before (either S1, or S2 in case of S1S2S3)
        has already generated all the LATs that it is going to before mintime, otherwise previous beats might get counted after mintime.
        Is this actually working properly, then? Even for S1S2S3 protocol?
        It is probably safe to assume that as long as the upstroke has crossed after minTime, that all LATs after mintime are the new beat.
        This might not be true if the wave hasn't crossed the strip in totality yet... could check for monotonous increase?

    """

    PRINT = True

    CV, APD = np.nan, np.nan

    minTime = S1 * (num_stim-1) + S2  # NOTE: this allows for using the S1 mode to calculate CV(S2), as long as S2 is given (defaults is zero)
    lat_list, thresh_list = ["act", "rep"], ["0.7", "0.1"]

    # check strip length, so we know where to measure CV and APD
    L = X[:,1].max() / 4.0 # represents an 8th of the distance
    CONST = 3.5 # used to exclude vertices very near stimulus when consider if all vertices have an upstroke
    GOOD_SIZE = (X[:,1] > (-CONST * L)).sum()
    #print("GOOD SIZE:", GOOD_SIZE)
    #print("SIZE:", X[:,1].size)
    #input()

    try:

        for lat, thresh in zip(lat_list, thresh_list):

            #{{{ load results files
            if mode == "S1": # S1 mode, for CV, APD
                datafile = job_dir + "/t" + lat + "_" + thresh + "-thresh.dat"
            else: # S2 mode, for ERP
                datafile = job_dir + "/t" + lat + "_" + str(S2) + "-thresh.dat"

            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=UserWarning, append=1)
                data = np.loadtxt(datafile)

            # if file is empty, then no activations occurred
            if data.size == 0:

                if PRINT: print("[WARNING]: Raising value error because empty data file, for S1 pacing")
                raise ValueError

            #}}}

            #{{{ keep activation times after the last stimulus, check how many activations occured

            df = pd.DataFrame(data, columns = ["vertex", "lat"])
            #print("mode:", mode)
            #print(df.groupby(['vertex']).size())
            df = df[df["lat"] > minTime]
            df.sort_values(["vertex", "lat"], inplace = True)
            #new_cond = X[df.vertex.values.astype(np.int32), 1] > -CONST * L
            #print(df.groupby(['vertex']).size())
            #print("df:", df)
            #input()

            #plt.scatter(df[new_cond].vertex.values, df[new_cond].lat.values, color = "blue")
            #plt.scatter(df[~new_cond].vertex.values, df[~new_cond].lat.values, color = "red")
            #plt.title(lat)
            #plt.show()

            if True:  df.drop_duplicates(subset = "vertex", keep = "last", inplace = True)
            else:
                #{{{ estimate of S1_local that used last and second to last activations
                LAT_second_to_last = \
                    (df.groupby('vertex', as_index=False).apply(lambda x: x if len(x)==1 else x.iloc[[-2]]).reset_index(level=0, drop=True))["lat"]

                df.drop_duplicates(subset = "vertex", keep = "last", inplace = True)

                diff = df.lat.values - LAT_second_to_last
                plt.title("LAT difference between last and second-to-last activation ({:s})".format(lat))
                plt.xlabel("X")
                plt.ylabel("LAT diff")
                plt.scatter(X[:,1], diff)
                plt.show()
                #}}}

            #}}}

            #{{{ check whether all activations occurred
            #if df.shape[0] < X.shape[0]: # considering al activations (just slighly concerned that at stimulus points, upstroke may not get recorded)
            if df[ X[df.vertex.values.astype(np.int32),1] > -3.5 * L ].shape[0] < GOOD_SIZE: # just considering activations beyond small region near stimulus

                if PRINT: print("[WARNING]: Only {:d}/{:d} activations occurred for S1, setting NaN results.".format(df.shape[0], X.shape[0]))
                raise ValueError

            # get the LATs for the interpolation resolution mesh
            if lat == "act": LAT_dep = df.lat.values - minTime
            else:            LAT_rep = df.lat.values - minTime

            #}}}

        #{{{ for S1, do CV and APD calculations
        if True:

            # rotate the geometry for square simulation
            if geometry == "square":
                Rot = Rotation.from_euler('z', 45, degrees = True)
                X = Rot.apply(X)

            #{{{ plot the CV and APD on the entire mesh
            if False:
                import heartplot as hp
                if lat == "rep":
                    hp.plot_scalars_on_mesh(X, Tri, scalars = LAT_rep - LAT_dep, reverse = False)
            #}}}

            # get vertices where we had recordings
            vertices = df.vertex.values.astype(np.int32)

            # for SQUARE, get vertices only by the diagonal for 'square'
            if geometry == "square": vertices = np.argwhere((X[:,0] > -5.0) & (X[:,0] < 5.0)).flatten()

            #{{{ plot CV along the strip to see how much it varies with distance from stimulus and ends of the strip
            if False and geometry == "strip":
                # NOTE: resolution is 100, we'll use spacing of 200 just in case

                miny = X[:,1].min()
                maxy = X[:,1].max()

                CV_list = []
                for i in np.arange(miny, maxy, 300):
                    LAT_DEP_1 = LAT_dep[X[vertices,1] == i ]
                    #if PRINT: print(LAT_DEP_1)
                    LAT_DEP_1 = LAT_DEP_1[LAT_DEP_1 > 0].mean()

                    LAT_DEP_2 = LAT_dep[X[vertices,1] == (i + 300) ]
                    #if PRINT: print(LAT_DEP_2)
                    LAT_DEP_2 = LAT_DEP_2[LAT_DEP_2 > 0].mean()

                    gradLAT = (LAT_DEP_2 - LAT_DEP_1) / 300.0
                    CV = 1.0 / gradLAT
                    CV = CV/1000.0
                    CV_list.append(CV)
                    #if PRINT: print("Position: {:f}".format((i + 150)))
                    #if PRINT: print("CV: {:f} m/s".format(CV))

                plt.scatter(np.arange(miny, maxy, 300), CV_list)
                plt.ylim(0,30)
                plt.xlim(miny,maxy)
                plt.xlabel("position")
                plt.ylabel("CV")
                plt.show()
            #}}}

            #{{{ calculate CV & APD by fitting LAT to position
            if True:

                condition = (X[vertices, 1] > -L) & (X[vertices, 1] < +L)

                xval = X[vertices,1][condition]
                #if PRINT: print("xval:", xval)

                # CV calculation
                # --------------
                LAT_ups = LAT_dep[vertices][condition]

                # fit LAT = m*x + c to get m
                A = np.vstack([np.ones(xval.shape[0]), xval]).T
                c, m = np.linalg.lstsq(A, LAT_ups, rcond = None)[0]
                LAT_pred = c + m*xval
                #plt.scatter(xval, LAT_ups, zorder = 0); plt.plot(xval, LAT_pred, zorder = 1, color = "red", linestyle = "--");
                #plt.title("LAT"); plt.xlabel("X"); plt.show() # plot LAT along the strip
                CV = (1.0/m) / 1000.0
                #if PRINT: print("CV: {:f} m/s".format(CV))

                # APD calculation
                # ---------------
                LAT_downs = LAT_rep[vertices][condition]

                APD = (LAT_downs - LAT_ups)
                APD = APD.mean()

                #if PRINT: print("APD: {:f} ms".format(APD))
            #}}}

        #}}}

    except ValueError as e: CV, APD = np.nan, np.nan

    except IndexError as e: CV, APD = np.nan, np.nan

    if mode == "S1":
        print("results (S1: {:d}) CV, APD :: {:1.3f} {:3.1f}".format(S1, CV, APD))
    if mode == "S2":
        print("results (S1: {:d}, S2: {:d}) CV, APD :: {:1.3f} {:3.1f}".format(S1, S2, CV, APD))

    return CV, APD

#}}}


def logistic_func(x, a, b, c):
    """Logistic function."""
    return a * (1 - b*np.exp( -x/c ) )


#{{{ parser
def parser():
    # Generate the standard command line parser
    parser = tools.standard_parser()

    # force option "--overwrite-behaviour" (dest: 'overwrite_behaviour') to have default 'delete'
    for action in parser._actions:
        if action.dest == "overwrite_behaviour":
            action.default = "delete"


    group  = parser.add_argument_group('experiment specific options')
    # Add arguments
    group.add_argument('--nbeats',
                        type=int,
                        default=8,
                        help='Number of beats for S1 pacing at CI1')
    group.add_argument('--mesh',
                        default='strip',
                        choices=['strip',
                                 'square'],
                        help='Which geometry to use.')
    group.add_argument('--S1',
                        type=int,
                        default=240,
                        choices=[600, 550, 500, 450, 400, 350, 300, 240],
                        help='S1 pacing value.')
    group.add_argument('--S2',
                        type=int,
                        default=0,
                        help='S2 pacing value for S1, ..., S1, S2, S3 protocol')
    group.add_argument('--filename',
                        type=str,
                        required = True,
                        help='Name of parameters file.')

    return parser
#}}}


def jobID(args):
    """
    Generate name of top level output directory.
    """
    today = date.today()
    return '{}_basic'.format(today.isoformat())


# supply extra arguments to run() via a keyword, which can be modified in __main__
def wrapper(fn):
    @functools.wraps(fn)
    def wrapped(*args, **kwargs):
        #try:
        return fn(*args, **kwargs, params = PARAMS)
        #except:
        #    print("[ERROR]: 'PARAMS' must be set in order to call run() function.")
        #    return None
    return wrapped


@tools.carpexample(parser, jobID, clean_pattern='^(\d{4}-\d{2}-\d{2})|(.txt)|(.dat)') #, simple = True)
@wrapper
def run(args, job, params = []):
    """ Run tissue strip simulation with given parameters.
    """

    print("params:", params)

    # construct basis carp command
    cmd = tools.carp_cmd('settings.par')

    # read mesh
    mesh = args.mesh + "/" + args.mesh
    X = np.loadtxt(mesh + ".pts", skiprows = 1)
    Tri = np.genfromtxt(mesh + ".elem", skip_header = 1, usecols = (1,2,3), dtype = np.int32)


    #{{{ construct the stimuli

    STRIP_START_X = X[:,0].min()
    STRIP_WIDTH = 2.0 * np.abs(STRIP_START_X)
    STRIP_START_Y = X[:,1].min()

    #{{{ S1 stimuli
    stim = [ '-num_stim',                   1,
             '-stimulus[0].name',           'S1',
             '-stimulus[0].stimtype',       0,
             '-stimulus[0].strength',       1.0,
             '-stimulus[0].duration',       1.0,
             '-stimulus[0].npls',           args.nbeats,
             '-stimulus[0].bcl',            args.S1,
             '-stimulus[0].x0',             STRIP_START_X,
             '-stimulus[0].xd',             STRIP_WIDTH,
             '-stimulus[0].y0',             STRIP_START_Y,
             '-stimulus[0].yd',             500.,
             '-stimulus[0].dump_vtx_file',  0
             ]
    #}}}

    #{{{ S2 stim for S3 protocol
    if args.S2 > 0:

        stim_for_S2 = \
                ['-num_stim',                   2, # this should overwrite num_stim = 1
                 '-stimulus[1].name',           'S2forS3',
                 '-stimulus[1].stimtype',       0,
                 '-stimulus[1].strength',       1.0,
                 '-stimulus[1].start',          (args.nbeats - 1) * args.S1 + args.S2,
                 '-stimulus[1].duration',       1.0,
                 '-stimulus[1].npls',           1,
                 '-stimulus[1].bcl',            1,
                 '-stimulus[1].x0',             STRIP_START_X,
                 '-stimulus[1].xd',             STRIP_WIDTH,
                 '-stimulus[1].y0',             STRIP_START_Y,
                 '-stimulus[1].yd',             500.,
                 '-stimulus[1].dump_vtx_file',  0
                ]

        stim = stim + stim_for_S2
    #}}}

    #}}}


    #{{{ set cell and conductivity parameters
    D, tau_in, tau_out, tau_open, tau_close = params
    G = 2.0 * D * 1000.
    mMS_par = 'V_gate=0.1,a_crit=0.1,tau_in={:f},tau_out={:f},tau_open={:f},tau_close={:f}'.format(tau_in, tau_out, tau_open, tau_close)

    param = [ '-imp_region[0].im_param',  mMS_par,
              '-gregion[0].g_el',         G,
              '-gregion[0].g_et',         G,
              '-gregion[0].g_il',         G,
              '-gregion[0].g_it',         G
            ]
    #}}}


    # S2(S3) values
    S2min, S2max = 50, 500


    #{{{ run S1 pacing simulation

    # save S1 state
    TSAVE = "{:.1f}".format((args.nbeats - 1) * args.S1 + args.S2 + S2min) # NOTE: including S2, for S1S2S3 pacing
    save = ['-num_tsav',      1,
            '-tsav[0]',       TSAVE,
            #'-write_statef',  job.ID + "/statefile" # NOTE: this does not work properly with carputils, it ignores requested name
           ]


    # construct full carp command
    cmd_S1 = cmd + stim + param + save

    # NOTE: probably require less time, but this is to be safe
    cmd_S1 += [ '-simID',     job.ID,
                '-meshname',  mesh,
                '-tend',      args.nbeats * args.S1 + args.S2,
                '-spacedt',   args.nbeats * args.S1 + args.S2,
                '-timedt',    args.nbeats * args.S1 + args.S2
              ]

    print("cmd:", cmd_S1)

    job.carp(cmd_S1)
    #job.mpi(cmd)

    # calculate results
    CV, APD = calc_CV_APD(job.ID, X, Tri, args.mesh, args.nbeats, S1 = args.S1, mode = "S1", S2 = args.S2) # NOTE: include S2 here

    #}}}


    # S2 restitution and ERP
    # ----------------------

    S2list = [] # list of tested values
    CVlist = [] # list of tested values
    APDlist = [] # list of tested values
    ERP = np.nan

    if ~np.isnan(CV) and ~np.isnan(APD):

        # S2 protocol, loading the saved state from the S1 protocol

        # S2low & S2high define search brackets (initially S2min & S2max)
        S2low, S2high, = S2min, S2max

        test_S2max = True # required to stop search if maxS2 is ERP
        while True:

            print("still in while")

            if test_S2max: S2 = S2high
            else:          S2 = int( (S2low + S2high)/2.0 )

            S2list.append(S2)
            #print("S2list", S2list)
            #input()

            #{{{ S2 stimuli
            stim_S2 = [ '-num_stim',                   1,
                        '-stimulus[0].name',           'S2',
                        '-stimulus[0].start',          (args.nbeats - 1) * args.S1 + args.S2 + S2,
                        '-stimulus[0].stimtype',       0,
                        '-stimulus[0].strength',       2.0,
                        '-stimulus[0].duration',       1.0,
                        '-stimulus[0].npls',           1,
                        '-stimulus[0].bcl',            1,
                        '-stimulus[0].x0',             STRIP_START_X,
                        '-stimulus[0].xd',             STRIP_WIDTH,
                        '-stimulus[0].y0',             STRIP_START_Y,
                        '-stimulus[0].yd',             500.,
                        '-stimulus[0].dump_vtx_file',  0
                      ]
            #}}}

            #{{{ construct full carp command for S2
            cmd_S2 = cmd + stim_S2 + param

            cmd_S2 += [ '-simID',         job.ID,
                        '-meshname',      mesh,
                        '-tend',          args.nbeats * args.S1 + args.S2 + S2,
                        '-spacedt',       args.nbeats * args.S1 + args.S2 + S2,
                        '-start_statef',  job.ID + "/state." + TSAVE,   # NOTE: works only because forced TSAVE to string with 1 decimal place
                        '-timedt',        args.nbeats * args.S1 + args.S2 + S2,
                        # NOTE: it will complain about not being able to read these files, probably because this continues from the S1 simulation
                        '-lats[0].ID',    "tact_" + str(S2),
                        '-lats[1].ID',    "trep_" + str(S2)
                        #'-lats[0].all',   0,
                        #'-lats[1].all',   0,
                      ]

            job.carp(cmd_S2)

            #}}}

            CV_S2, APD_S2 = calc_CV_APD(job.ID, X, Tri, args.mesh, args.nbeats, S1 = args.S1, mode = "S2", S2 = S2)

            CVlist.append(CV_S2)
            APDlist.append(APD_S2)

            # test if this value of S2 did not activate
            if np.isnan(CV_S2):
                ERP = S2


            #{{{ whether to continue searching for ERP

            # if testing the max S2 and no activation, break loop
            if test_S2max:
                if ~np.isnan(ERP): break
                test_S2max = False

            if np.isnan(ERP): # if this S2 was above ERP
                S2high = S2
            else: # if this S2 was below ERP
                S2low = S2

            if (S2high - S2low == 1) or (S2high - S2low == 0): # break loop because found ERP to 1ms resolution
                print("S2low, S2high:", S2low, S2high)
                ERP = S2low
                break
            else:
                ERP = np.nan # reset ERP back to nan, now that search range has been changed, and continue to refine search

            #}}}


    print("tested S2 values:", S2list)
    print("S2: CV values :", CVlist)
    print("S2: APD values:", APDlist)

    print("Finished!")
    print("CV, APD, ERP: {:1.3f}, {:3.1f}, {:3.1f}".format(CV, APD, ERP))

    S2array = np.array([S2list, CVlist, APDlist]).T
    print("S2array:\n", S2array)
    S2array = S2array[~np.isnan(S2array).any(1)]
    print("S2array:\n", S2array)

    # FIXME: the return should be CV, APD, ERP, and an array of S2
    return CV, APD, ERP, S2array


#{{{ function to append to results file
def save_results(results_file, CV, APD, ERP, CVparams, APDparams):

    with open(results_file, "a") as myfile:
        print("Appending results to {:s}".format(results_file))

        #myfile.write("{:1.3f} {:3.1f} {:3.1f}\n".format(CV, APD, ERP)) # no CV(S2) information

        myfile.write("{:1.3f} {:3.1f} {:1.1f} {:1.3f} {:1.3f} {:1.3f} {:1.3f} {:1.3f} {:1.3f} \n".format( \
                     CV, APD, ERP, \
                     CVparams[0],  CVparams[1],  CVparams[2], \
                     APDparams[0], APDparams[1], APDparams[2] ))
#}}}


if __name__ == '__main__':

    # In order to use the 'carpexample' framework, I have wrapped run() with a wrapper function that calls the original run() with argument param = PARAM.
    # The value PARAM must therefore be set for this to work.
    # For simplicity, I have added the option --filename to the parser and accessed this value here, in order to load the parameters file.


    #{{{ set parameters file using --filename argument
    prs = parser()
    args = prs.parse_args(sys.argv[1:])

    try:
        data = np.loadtxt(args.filename)
        print("Loaded parameters from file {:s}".format(args.filename))
    except OSError as e:
        print("[ERROR]: file {:s} not found".format(args.filename))
        exit()
    #}}}


    #{{{ name results file automatically
    results_file = os.path.split(args.filename)[0] + "/" \
                 + args.filename.split("/")[-1].split(".")[0] \
                 + "_S1" + str(args.S1) \
                 + "_S2" + str(args.S2) \
                 + ".dat"
                 #+ datetime.strftime(datetime.today(), '%a%d%b') \

    print("results filename:", results_file)

    # note: can access S1 and S2 from the filename using int(tmp.split("S1")[-1].split("_")[0].split(".")[0])
    print("S1:", int(results_file.split("S1")[-1].split("_")[0].split(".")[0]))
    print("S2:", int(results_file.split("S2")[-1].split("_")[0].split(".")[0]))
    #}}}


    # for recording all restitution curves
    S2_rest_all = []

    # loop over parameters
    for PARAMS in data:

        # run simulator and obtain results
        CV, APD, ERP, S2rest = run()

        try:

            if np.isnan(CV) or np.isnan(APD):

                CVparams   = np.full(3, np.nan)
                APDparams  = np.full(3, np.nan)

            else:

                print("CV, APD, ERP:", CV, APD, ERP)
                print("S2rest:", S2rest)

                # sort S2rest by S2 value
                S2rest = S2rest[ S2rest[:,0].argsort() ] # sort by S2
                S2_rest_all.append(S2rest) # append to S2_rest_all

                # discard inflection of CV(S2) at low S2
                REST = S2rest[np.argmin(S2rest[:,1]):, :]
                print("REST:", REST)

                # fit CV(S2) & APD(S2) restitutions with logistic function
                for TYPE in [1, 2]:

                    if TYPE == 1: bounds = ( [0, 0, 0], [1000, 1000, 1000] )
                    if TYPE == 2: bounds = ( [0, 0, 0], [1000, 1000, 1000] )


                    # fit to logistic
                    #params, cov = curve_fit(logistic_func, REST[:,0], REST[:,TYPE], bounds = ( [0,0,0,0], [1000, 1000, 1000, 1000] ))
                    params, cov = curve_fit(logistic_func, REST[:,0], REST[:,TYPE], bounds = bounds )
                    print("params:", params)

                    if TYPE == 1: CVparams  = params
                    if TYPE == 2: APDparams = params

                    if False: # plots

                        plt.scatter(REST[:,0], REST[:,TYPE])
                        S2 = np.linspace(REST[:,0].min(), REST[:,0].max(), 100, endpoint = True)
                        plt.plot(S2, logistic_func(S2, *params), color = "red", linestyle = "--")

                        if TYPE == 1:
                            TITLE = "CV restitution"
                        else:
                            TITLE = "APD restitution"

                        plt.title(TITLE)
                        plt.show()


            # call results saving routine to append line to results file
            save_results(results_file, CV, APD, ERP, CVparams, APDparams)

        except ValueError as e:
            print("[ERROR]:", e)
            save_results(results_file, np.nan, np.nan, np.nan, np.full(3, np.nan),  np.full(3, np.nan))

        except RuntimeError as e:
            print("[ERROR]:", e)
            save_results(results_file, np.nan, np.nan, np.nan, np.full(3, np.nan),  np.full(3, np.nan))


    # pickle the raw restitutions results to a temporary file
    tmpfile = '/tmp/restitutions.pkl'
    print("Saving resitutions to", tmpfile)
    with open(tmpfile, 'wb') as f:
        pickle.dump(S2_rest_all, f)
