# coding: utf-8
import numpy as np
import pandas as pd
import scipy as sp
import scipy.stats
import emcee
import random
from math import *

import subprocess
import pickle as pkl
import os
import shutil
from tqdm import tqdm

#import matplotlib
#import matplotlib.pyplot as plt
#import corner

from ratios import r02, r10, get_ratios

np.set_printoptions(precision=10)
np.random.seed(42) # for reproducibility 

Z_X_solar = 0.02307 # GS98

nwalkers, nthreads, niter, nperturb = 40, 8, 1000, 100000

# Priors and starting point 
# Normal on t_0 and reaction rates; flat on all others 
names     = [   't',   'M',  'Y',    'Z',     'a',    'b',  't_0', 'r1', 'r2', 'r3', 'r4', 'r5']
X_max     = [13.799,     1,  0.4,   0.02,     2.5,    0.2, 14.219,    2,    2,    2,    2,    2]
X_min     = [     0,   0.5,  0.2,  0.001,     0.2,   -0.2, 13.379,    0,    0,    0,    0,    0]
prior_mu  =                                               [13.799,    1,    1,    1,    1,    1]
prior_cov =                              np.diag(np.array([ 0.021, 0.01, 0.05, 0.05, 0.08, 0.07])**2)
start     = [  10.7, 0.725, 0.251,  0.0059,   1.90,    0., 13.799,   1.,   1.,   1.,   1.,   1.]
start_cov = [   1.2, 0.044, 0.033,  0.0013,    0.1,   0.01] + list(np.diag(prior_cov))
ndim = len(names)

# Directories and save files 
MCMC_dir = "MCMC_calcs"
save_dir = "MCMC_save"

if not os.path.exists(MCMC_dir):
    os.mkdir(MCMC_dir)

if not os.path.exists(save_dir):
    os.mkdir(save_dir)

freqs_filename       = "7970740-freqs.dat"
Sigma_filename       = os.path.join(save_dir, "Sigma.npy")
Sigma_inv_filename   = os.path.join(save_dir, "Sigma_inv.npy")
norm_factor_filename = os.path.join(save_dir, "norm_factor.npy")
state_filename       = os.path.join(save_dir, "state.pkl")
chain_filename       = os.path.join(save_dir, "chain.dat")

# Load stellar data 
Teff_solar   = 5772.
Teff_stellar = 5309.
e_Teff_stellar = 77.
FeH_stellar  = -0.54
e_FeH_stellar = 0.1

freqs = pd.read_table(freqs_filename, sep='\s+')

obs_r02 = r02(freqs)
obs_r10 = r10(freqs)

r02_interp_freqs = obs_r02['freqs']
r10_interp_freqs = obs_r10['freqs']

print('r02', r02_interp_freqs)
print('r10', r10_interp_freqs)

y_star = pd.DataFrame([[Teff_stellar/Teff_solar, FeH_stellar] + \
        list(obs_r02['ratios']) + list(obs_r10['ratios'])],
    columns=['Teff', 'Fe_H'] + list(obs_r02['names'] + obs_r10['names']))

# log likelihood calculation: 
# ln L = - [ ln(|Sigma|) + k*ln(2*pi) + chi2 ] / 2. 
# where |Sigma| is the determinant of the covariance matrix 
# and k is the dimensionality of the input data 

# Now perform 100,000 Monte Carlo realizations to obtain Sigma 
if os.path.isfile(Sigma_inv_filename) and os.path.isfile(norm_factor_filename):
    Sigma_inv = np.load(Sigma_inv_filename)
    norm_factor = np.load(norm_factor_filename)
else:
    print('Computing ratios')
    r02s = get_ratios(r02, freqs, n_perturb=nperturb, seed=0)
    r10s = get_ratios(r10, freqs, n_perturb=nperturb, seed=0)
    rs = np.hstack((r02s, r10s))
    
    Sigma = np.cov(rs.T)
    Sigma = np.vstack([np.zeros(Sigma.shape[0]), 
                       np.zeros(Sigma.shape[0]), Sigma])
    Sigma = np.insert(Sigma, 0, np.zeros(Sigma.shape[0]), axis=1)
    Sigma = np.insert(Sigma, 0, np.zeros(Sigma.shape[0]), axis=1)
    Sigma[0,0] = (e_Teff_stellar/Teff_solar)**2
    Sigma[1,1] = e_FeH_stellar**2
    print('Condition number:', np.linalg.cond(Sigma))
    
    Sigma_inv = np.linalg.inv(Sigma)
    norm_factor = np.log(np.linalg.det(Sigma)) + len(y_star) * np.log(2*np.pi)
    np.save(Sigma_filename,       Sigma)
    np.save(Sigma_inv_filename,   Sigma_inv)
    np.save(norm_factor_filename, norm_factor)

def lnprior(theta, X_max=X_max):
    _theta = np.copy(theta)
    X_max[names.index('t')] = _theta[names.index('t_0')]
    if any(_theta >= X_max) or any(_theta <= X_min):
        print("Out of bounds")
        return -np.inf
    #return 0.0 # flat prior 
    _lnprior = sp.stats.multivariate_normal.logpdf(
        _theta[-len(prior_mu):], mean=prior_mu, cov=prior_cov)
    return _lnprior 

def lnlike(theta):
    _theta = np.copy(theta)
    
    age = _theta[names.index('t')]
    
    _theta[names.index('t')]   = float(str(_theta[names.index('t')])   + 'e9')
    _theta[names.index('t_0')] = float(str(_theta[names.index('t_0')]) + 'e9')
    
    # generate a random temp filename and make sure it's not in use 
    while True:
        rand_bits = str(random.getrandbits(64))
        tmp_dir = os.path.join(MCMC_dir, rand_bits)
        if not os.path.exists(tmp_dir):
            break
    
    # call ASTEC 
    bash_cmd = "./astec_mcmc.sh -r -d " + MCMC_dir + " -n " + rand_bits + \
        ' -'.join([''] + list(map(lambda x,y: x+' '+str(y), names, _theta)))
    print(bash_cmd)
    process = subprocess.Popen(bash_cmd.split(), shell=False)
    try: 
        process.wait(timeout=6 * 60 * 60) # seconds = 6 hours 
    except subprocess.TimeoutExpired:
        process.terminate()
    
    freq_file = tmp_dir + '.freqs'
    hist_file = tmp_dir + '.dat'
    gong_file = tmp_dir + '.FGONG.dat'
    if not os.path.exists(freq_file) or \
       not os.path.exists(hist_file) or \
       not os.path.exists(gong_file):
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir, ignore_errors=True)
        print(rand_bits, 'Files not generated')
        return -np.inf
    
    nus = pd.read_table(freq_file, sep='\s+', comment='#', 
        names=['l', 'n', 'nu', 'E', 'n_p', 'n_g'])
    
    # compute ratios from frequencies 
    r02s = r02(nus, interp_freqs=r02_interp_freqs)
    r10s = r10(nus, interp_freqs=r10_interp_freqs)
    
    # obtain Teff from history file 
    hist_file = tmp_dir + '.dat'
    hist = pd.read_table(hist_file, sep='\s+', comment='#', 
        names=['Model', 'M', 'age', 'R', 'Teff', 'L', 'X_c', 'qc'])
    Teff = hist['Teff'].values[-1]
    
    # obtain [Fe/H] from FGONG file 
    gong_file = tmp_dir + '.FGONG.dat'
    gong = pd.read_table(gong_file, sep='\s+')
    FeH = log10(gong['Z'].values[0] / gong['X'].values[0] / Z_X_solar)
    
    # delete files 
    os.remove(freq_file)
    os.remove(hist_file)
    os.remove(gong_file)
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir, ignore_errors=True)
    
    # check that the age is close to the requested age 
    if np.abs(hist['age'].values[-1] - age)/age > 0.1:
        print(rand_bits, 'Age mismatch')
        return -np.inf 
    
    # check that we were able to compute the ratios 
    if r02s['ratios'] is None or r10s['ratios'] is None:
        print(rand_bits, 'Ratios were not computed')
        return -np.inf
    
    for name in obs_r02['names']:
        if name not in r02s['names']:
            print(rand_bits, name, 'not computed')
            return -np.inf
    
    for name in obs_r10['names']:
        if name not in r10s['names']:
            print(rand_bits, name, 'not computed')
            return -np.inf
    
    # compute likelihood 
    y = pd.DataFrame([[Teff/Teff_solar, FeH] + \
            list(r02s['ratios']) + list(r10s['ratios'])],
        columns=['Teff', 'Fe_H'] + list(obs_r02['names'] + obs_r10['names']))
    
    if y.isnull().values.any():
        print(rand_bits, "NaNs")
        return -np.inf
    
    resid = y.iloc[0] - y_star.iloc[0]
    chi2 = np.dot(resid, np.dot(Sigma_inv, resid))
    return -(norm_factor + chi2)/2.

def lnprob(theta, i=0, N=0, size=0):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta)

if os.path.isfile(state_filename): 
    # if the calculation has already started, continue from there 
    print("Resuming calculation")
    state_file = open(state_filename, 'rb')
    pos, prob, state = pkl.load(state_file)
    state_file.close()
else: 
    print("Initializing calculation")
    # otherwise, start in a small ball around a good initial guess 
    prob, state = None, None
    
    pos = np.array([[np.random.normal(start[ii], start_cov[ii]/10) 
                     for ii in range(ndim)]
           for jj in range(nwalkers)])
    
    chain_file = open(chain_filename, 'w')
    chain_file.close()

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, threads=nthreads)
#chain = np.empty((nwalkers, niter, ndim))
#flatchain = np.empty((0, ndim))

for ii in tqdm(range(niter)): 
    for results in sampler.sample(pos, lnprob0=prob, rstate0=state, 
            storechain=False):
        pass
    pos, prob, state = results 
    
    chain_file = open(chain_filename, 'a')
    for k in range(pos.shape[0]):
        chain_file.write("{0:4d} {1:s}\n".format(k, " ".join(map(str, pos[k]))))
    chain_file.close()
    
    state_file = open(state_filename, 'wb')
    pkl.dump(results, state_file)
    state_file.close()
    
    #flatchain = np.concatenate((flatchain, pos))
    #chain[:, ii, :] = pos 


#chain_pkl = open('chain.pkl', 'wb')
#pkl.dump(sampler.chain, chain_pkl)
#chain_pkl.close()
