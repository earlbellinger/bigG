import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
from corner import corner
import numpy as np

DF = pd.read_table('MCMC_save/chain.dat', sep='\s+', header=None,
    names=['index', 'age', 'M', 'Y', 'Z', 'alpha', 'beta', 't0', 
           'rate1', 'rate2', 'rate3', 'rate4', 'rate5'])

DF = DF.iloc[-100000:,]
print("shape:", DF.shape)
#corner(DF.iloc[-70000:,1:])
corner(DF.iloc[:,1:])
plt.savefig('bigG_corner.pdf')



print("beta")
print(DF['beta'].mean())
print(DF['beta'].std())
print()
print()
print("t0")
print(DF['t0'].mean())
print(DF['t0'].std())
print()
print()
print("dot G")
print((DF['beta']/(DF['t0']*10**9)).mean())
print((DF['beta']/(DF['t0']*10**9)).std())

vals = [4.01, 5.21, 0.56, 2.08, 1.66]
DF.iloc[:,-5:].mean() #* vals 
DF.iloc[:,-5:].std() #* vals 


DF.iloc[:,1:6].mean()
DF.iloc[:,1:6].std()


X_max     = [14.219,    2,    2,    2,    2,    2]
X_min     = [13.379,    0,    0,    0,    0,    0]
prior_mu  =         [13.799,    1,    1,    1,    1,    1]
prior_cov = np.diag([ 0.021, 0.01, 0.05, 0.05, 0.08, 0.07])

X_labels = ['t_0', 'Rate 1', 'Rate 2', 'Rate 3', 'Rate 4', 'Rate 5']
ndim = len(X_labels)
arr = DF.iloc[:,-6:]
fig = corner(arr, labels=X_labels,
    #quantiles=[0.5],
    show_titles=True,
    title_fmt='.4f',
    #color='b',
    truths=np.percentile(arr, 50, axis=0),
    #truth_color='b',
    title_kwargs={"fontsize": 12})

axes = np.array(fig.axes).reshape((ndim, ndim))

for ii in range(len(X_labels)):
    ax = axes[ii, ii]
    xs = np.linspace(X_min[ii], X_max[ii], num=100)
    ys = sp.stats.norm.pdf(xs, loc=prior_mu[ii], scale=prior_cov[ii, ii])
    ymax = ax.get_ylim()[1]
    ax.plot(xs, ymax * ys/max(ys) * 0.92, 'r--')

fig.savefig('Sfactor-corner.pdf')
fig.show()