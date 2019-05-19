#!usr/bin/env python
#### Call dispatch.sh with quasi-random inputs 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

import numpy as np
import subprocess
import argparse
from time import sleep

from sys import path
path.append('/home/earl/asteroseismology/scripts')
from sobol_lib import i4_sobol

def main(arguments):
    parser = argparse.ArgumentParser()
    init = parser.add_argument_group('initial conditions')
    init.add_argument('-M', default=[0.6, 0.9], nargs=2, type=float,
                      help='range of masses (solar = 1)')
    init.add_argument('-Y', default=[0.22, 0.30], nargs=2, type=float, 
                      help='range of helium values (solar = 0.273)')
    init.add_argument('-Z', default=[0.001, 0.0185], nargs=2, type=float,
                      help='range of metallicity values (solar = 0.0185)')
    init.add_argument('-a', '--alpha', default=[0.5, 2.5], nargs=2, type=float,
                      help='range of mixing length parameter values (solar = 1.83)')
    init.add_argument('-t', '--age', default=[5, 13.799], nargs=2, type=float,
                      help='range of ages')
    init.add_argument('-b', '--beta', default=[-0.1, 0.1], nargs=2, type=float,
                      help='range of gravitational power law parameters')
    #init.add_argument('-l', '--logs', default=[0, 0, 0, 0, 0, 0], nargs='6',
    #                  type=list, 
    #                  help='booleans of whether to log M, Y, Z, alpha, age, beta')
    
    rates = parser.add_argument_group('rates')
    rates.add_argument('--rate1', default=[0.9, 1.1], nargs=2, type=float,
                      help='')
    rates.add_argument('--rate2', default=[0.5, 1.5], nargs=2, type=float,
                      help='')
    rates.add_argument('--rate3', default=[0.5, 1.5], nargs=2, type=float,
                      help='')
    rates.add_argument('--rate4', default=[0.2, 1.8], nargs=2, type=float,
                      help='')
    rates.add_argument('--rate5', default=[0.3, 1.7], nargs=2, type=float,
                      help='')
    
    job = parser.add_argument_group('job')
    job.add_argument('-N', '--tracks', default=65536, type=int, 
                     help='number of tracks to generate')
    job.add_argument('-s', '--skip', default=20000, type=int,
                     help='offset for sobol numbers')
    job.add_argument('-d', '--directory', default="simulations", type=str,
                     help='offset for sobol numbers')
    job.add_argument('-r', '--remove', default=False, action='store_true',
                     help='delete models upon completion')
    
    physics = parser.add_argument_group('physics')
    physics.add_argument('-c', '--chem_ev', default=False, action='store_true',
        help='set Y = c*Z + 0.2463, where the -Y flag becomes the c range (solar c = 1.4276221707417)')
    physics.add_argument('-A', '--abundances', default=3, type=int,
        help='change mixture of Z (3=GS98, 6=AGSS09)')
    
    args = parser.parse_args(arguments)
    #print(args)
    
    ranges = np.vstack((args.M, args.Y, args.Z, args.alpha, args.age, 
        args.beta, args.rate1, args.rate2, args.rate3, args.rate4, args.rate5))
    
    #for i in range(len(args.logs)):
    #    if args.logs[i]:
    #        ranges[i] = np.log10(ranges[i])
    
    dispatch(ranges=ranges, tracks=args.tracks, 
        #logs=args.logs, 
        directory=args.directory, remove=args.remove, 
        skip=args.skip, chem_ev=args.chem_ev, abundances=args.abundances)

def dispatch(ranges, tracks, #logs, 
        directory, remove=0, skip=0, chem_ev=0,
        abundances=3):
    shift = ranges[:,0]
    scale = np.array([(b-a) for a,b in ranges])
    init_conds = []
    for i in range(skip, tracks+skip):
        vals = shift+np.array(i4_sobol(len(ranges), i)[0])*scale
        #for j, val in enumerate(vals):
        #    if logs[j]:
        #        vals[j] = 10**val if not np.isnan(val) else 0
        init_conds += [[tmp for tmp in vals]]
        for j, val in enumerate(vals):
            if np.isnan(vals[j]):
                vals[j] = 0
        if chem_ev:
            vals[1] = vals[1] * vals[2] + 0.2463
        
        bash_cmd = "./astec_dispatch.sh -d %s "\
            "-n %d "\
            "-M %.8f -Y %.8f -Z %.8f -a %.8f -t %.8fe9 -b %.8f "\
            "-r1 %.8f -r2 %.8f -r3 %.8f -r4 %.8f -r5 %.8f "\
            "-A %d "\
            "%s"%\
            tuple([directory] + 
                  [i] + 
                  [val for val in vals] + 
                  [abundances] +
                  ["-r " if remove else ""])
        print(bash_cmd)
        #exit()
        process = subprocess.Popen(bash_cmd.split(), shell=False)
        process.wait()
    #np.savetxt('initial_conditions.dat', np.array(init_conds))

if __name__ == '__main__':
    import sys
    exit(main(sys.argv[1:]))

