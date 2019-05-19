import numpy as np 
from scipy.interpolate import interp1d

def r02(nus, interp_freqs=None, perturb=False, verbose=False):
    """r_02(n) = (nu_{n,0} - nu_{n-1, 2}) / 
                 (nu_{n,1} - nu_{n-1, 1})
       where nus is a pandas DataFrame with columns l, n, nu (, dnu)"""
    
    nus_ = nus.copy()
    
    if perturb:
        nus_['nu'] = np.random.normal(nus_['nu'], nus_['dnu'])
    
    ell0, ell1, ell2 = (nus_[nus_['l'] == ell] for ell in range(3))
    
    ratios = []
    freqs  = []
    names  = []
    modes  = []
    for idx, mode in ell0.iterrows():
        n = mode['n'] # n, 0
        
        ell1m1 = ell1['n'] == (n - 1) # n-1, 1
        ell1n0 = ell1['n'] ==  n      # n,   1
        ell2m1 = ell2['n'] == (n - 1) # n-1, 2
        ell2n0 = ell2['n'] ==  n      # n,   2
        
        if not (sum(ell1m1) and sum(ell1n0) and sum(ell2m1) and sum(ell2n0)):
            if verbose:
                print('skipping', n)
            continue 
        
        num = float(mode['nu'])         - float(ell2[ell2m1]['nu']) # nu_{n,0} - nu_{n-1, 2}
        den = float(ell1[ell1n0]['nu']) - float(ell1[ell1m1]['nu']) # nu_{n,1} - nu_{n-1, 1}
        ratio = num/den
        
        ratios += [ratio]
        freqs  += [float(ell2[ell2n0]['nu'])] # nu_{n, 2}
        names  += ["r02_"  + str(int(n))]
        modes  += ["nu_0_" + str(int(n)), 
                   "nu_1_" + str(int(n-1)),
                   "nu_1_" + str(int(n)),
                   "nu_2_" + str(int(n-1)),
                   "nu_2_" + str(int(n))]
        
        if verbose:
            print('n:', n, 'num/den', num/den)
    
    if interp_freqs is None: 
        interp_freqs = freqs 
    
    if verbose:
        print(names)
        print(set(modes))
        print(interp_freqs)
    
    return {'names':  names,
            'modes':  set(modes),
            'freqs':  np.array(interp_freqs), 
            'ratios': interp1d(freqs, ratios, 
                              kind='cubic', 
                              bounds_error=False,
                              fill_value=np.nan)(interp_freqs)}

def r10(nus, interp_freqs=None, perturb=False, verbose=False):
    """r_10(n) = (     nu_{n-1, 1} 
                   - 4*nu_{n,   0} 
                   + 6*nu_{n,   1} 
                   - 4*nu_{n+1, 0}
                   +   nu_{n+1, 1}
              ) / -[8*(nu_{n+1, 0} 
                     - nu_{n,   0})]
       where nus is a pandas DataFrame with columns l, n, nu (, dnu)"""
    
    nus_ = nus.copy()
    
    if perturb:
        nus_['nu'] = np.random.normal(nus_['nu'], nus_['dnu'])
    
    ell0, ell1 = (nus_[nus_['l'] == ell] for ell in range(2))
    
    ratios = []
    freqs  = []
    names  = []
    modes  = []
    for idx, mode in ell0.iterrows():
        n = mode['n']
        
        ell1m1 = ell1['n'] == (n - 1)
        ell1n0 = ell1['n'] ==  n
        ell1p1 = ell1['n'] == (n + 1)
        
        ell0m1 = ell0['n'] == (n - 1)
        ell0p1 = ell0['n'] == (n + 1)
        
        if not (sum(ell1m1) and sum(ell1n0) and sum(ell1p1) \
                and sum(ell0m1) and sum(ell0p1)):
            if verbose:
                print('skipping', n)
            continue 
        
        num =     float(ell1[ell1m1]['nu']) \
              - 4*float(mode['nu'])         \
              + 6*float(ell1[ell1n0]['nu']) \
              - 4*float(ell0[ell0p1]['nu']) \
              +   float(ell1[ell1p1]['nu'])
        den = -8*(float(ell0[ell0p1]['nu']) \
                - float(mode['nu']))
        
        ratios += [num/den]
        freqs  += [float(mode['nu'])]
        names  += ["r10_" + str(int(n))]
        modes  += ["nu_1_" + str(int(n-1)), 
                   "nu_0_" + str(int(n)),
                   "nu_1_" + str(int(n)),
                   "nu_0_" + str(int(n+1)),
                   "nu_1_" + str(int(n+1))]
        
        if verbose:
            print('n:', n, 'num/den', num/den)
    
    if interp_freqs is None: 
        interp_freqs = freqs 
    
    return {'names':  names,
            'modes':  modes,
            'freqs':  np.array(interp_freqs), 
            'ratios': interp1d(freqs, ratios, 
                              kind='cubic', 
                              bounds_error=False,
                              fill_value=np.nan)(interp_freqs)}

def get_ratios(func, freqs, n_perturb=100000, seed=None):
    if seed is not None:
        np.random.seed(seed)
    return np.array((func(freqs, perturb=True)['ratios']
        for ii in range(n_perturb)))

#r02s = get_ratios(r02, freqs)
#r10s = get_ratios(r10, freqs)
