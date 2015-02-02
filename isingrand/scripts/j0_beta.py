#!/usr/bin/python
import numpy as np
from scipy.special import j0,jv
from scipy.optimize import fsolve
from matplotlib import rcParams
from matplotlib.pylab import *
from prettytable import PrettyTable

approx_roots = np.array([4.3,6.8,10.4,13.1,16.7,19.4])
x=np.linspace(0,20,100)

def beta(x):
    ints = range(1,101,2)
    return 2.0 * np.sum(jv(ints,x)/ints)

blist=[]
for xi in x:
    blist.append(beta(xi))
blist = np.array(blist)

rcParams.update({'font.size': 22})
xlabel(r'$\eta$')
plot(x,blist)
plot(x,j0(x))
plot(x,np.zeros(100))
legend([r'$\mathcal{J}_0(\eta)$',r'$\beta(\eta)$'],loc='upper right')
show()

roots=[]
for trial in approx_roots:
    roots.append(fsolve(beta,trial)[0])
roots=np.array(roots)    

tableofroots = PrettyTable(padding_width=2)
print
print "Roots of the beta function:"
print
tableofroots.add_column("Root #", np.arange(roots.size)+1)
tableofroots.align["Root #"] = "l" # Left align 
tableofroots.add_column("Root", roots)
tableofroots.add_column("J_0(Root)", j0(roots))
tableofroots.align["J_0(Root)"] = "r" # Left align 
print tableofroots
