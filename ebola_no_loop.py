# -*- coding: utf-8 -*-
#!/usr/bin/env python

#***********************************************************************
# Importing everything python needs in order to be smart.
#***********************************************************************
import numpy as np                  # Numpy tells python how to do array operations.
import matplotlib.pyplot as plt     # Matplotlib tells python how to make pretty plots.
import scipy.integrate as spi       # Python needs to know how to integrate systems of differential equations.
from scipy.integrate import odeint  # Teaches Python how to solve ODE's.

'''
Solving a modified SEIJR model (SEIJRDC) to simulate an Ebola virus outbreak.
'''

#***********************************************************************
# Define Functions
#***********************************************************************
#***********************************************************************
# Set of stiffly coupled differential equations.
# Solving the system dy/dt = f(y, t), where y = [S, E, I, J, R, D, C]
#***********************************************************************
def f(y, t):
    Si = y[0]
    Ei = y[1]
    Ii = y[2]
    Ji = y[3]
    Ri = y[4]
    Di = y[5]
    #Ci = y[6]
    dS = -beta   *  Si      * (Ii+l*Ji)     / N                      
    dE =  beta   *  Si      * (Ii+l*Ji)     / N       - k       * Ei
    dI =  k      *  Ei      - (alpha+gamma) * Ii                      
    dJ =  alpha  *  Ii      - gamma_r       * Ji                      
    dR =  gamma  *  (1.0-F) * Ii            + gamma_r * (1.0-F) * Ji
    dD =  gamma  *  F       * Ii            + gamma_r * F       * Ji
    #dC =  beta   *  Si      * (Ii+l*Ji)     / N          
    #return [dS, dE, dI, dJ, dR, dD, dC]
    return [dS, dE, dI, dJ, dR, dD]

#***********************************************************************
# End Define Functions
#***********************************************************************

#***********************************************************************
# Define constants - Set the disease parameters.
#***********************************************************************
# True constants.
#***********************************************************************
N       = 1000000.0             # Initial Population size.
gamma   = 1.0/5.0               # Average infectious period.
gamma_r = 1.0/7.0               # Average hospital stay.

#***********************************************************************
# Constants that need iterated.
#***********************************************************************
R_0     = 1.5                                          # Basic reprpoductive ratio of the virus.
F       = 0.5                                          # Case fatality proportion.

#***********************************************************************
# Set the initial conditions.
#***********************************************************************
E0 = 10.0                               # Initial number of exposed individuals.
I0 = 10.0                               # Initial number of infectious individuals in community.
J0 = 0.0                                # Initial number of hospitalized individuals.
R0 = 0.0                                # Initial number of individuals removed from isolation after recovery.
D0 = 0.0                                # Initial number of individuals removed from isolation after disease-induced death.
C0 = 0.0                                # Initial number of cumulative second cases.
S0 = N - (E0 + I0 + J0 + R0)            # Initial number of susceptible individuals.
#y0 = [S0, E0, I0, J0, R0, D0, C0]   # initial condition vector.
y0 = [S0, E0, I0, J0, R0, D0]   # initial condition vector.

#***********************************************************************
# Set the minimum and maximum time in [days] to simulate.
# Set the time grid. Start at day 0 and run until day 500 and do so 
# in 5000 steps (50 steps for each day).
#***********************************************************************
t  = np.linspace(0, 500, 50000)
#t  = np.linspace(0, 1000, 100000)

#***********************************************************************
# Effectiveness of isolation strategy (l=1 denotes no isolation; l=0 denotes perfect isolation).
# Determine l according to a given value of isol.
#***********************************************************************
l = 0.2

#***********************************************************************
# Average time from onset of symptoms to hispitalization.
# Determine alpha according to a given value of tdiag.
#***********************************************************************
alpha = 1.0/1.0

#***********************************************************************
# Latency period.
# Determine k according to a given value of incub.
#***********************************************************************
k = 1.0/8.0

#***********************************************************************
# Transmission Rate.
# Determine beta according to a given value of R0, alpha, gamma, gamma_r, and l.
#***********************************************************************
beta = (R_0 * (alpha + gamma) * gamma_r) / (alpha * l + gamma_r)
print beta
#***********************************************************************
# Solve the Differential Equations.
#***********************************************************************
soln = odeint(f, y0, t)
S    = soln[:, 0]      # Susceptible individuals.
E    = soln[:, 1]      # Exposed individuals.
I    = soln[:, 2]      # Infectious and symptomatic individuals.
J    = soln[:, 3]      # Hospitalized individuals.
R    = soln[:, 4]      # Individuals removed from isolation after recovery.
D    = soln[:, 5]      # Individuals removed from isolation after disease-induced death.
#C    = soln[:, 6]      # Cumulative second cases.

#***********************************************************************
# Plot the results.
#***********************************************************************
fig = plt.figure(figsize=[12.15, 7.06])
ax = fig.add_subplot(111)
ax.set_title('Ebola Virus Outbreak [Population vs. Time]' '\n' 'R0=' + str(R_0) + '  incub=' + str(k) + '  fat=' + str(F) + '  tdiag=' + str(alpha) + '  isol=' + str(l))
ax.set_ylabel('Population [# of People]')
ax.set_xlabel('Time Since Initial Outbreak [days]')
#ax.set_yticks(np.arange(0, 1000001, 100000))
#ax.semilogy(t, S, linewidth=2.0) 
#ax.semilogy(t, I, linewidth=2.0) 
#ax.semilogy(t, J, linewidth=2.0)
#ax.semilogy(t, R, linewidth=2.0)
#ax.semilogy(t, D, '--', linewidth=2.0) 
#ax.semilogy(t, C, linewidth=2.0)
ax.plot(t, S, linewidth=2.0)
ax.plot(t, E, linewidth=2.0)
ax.plot(t, I, linewidth=2.0)
ax.plot(t, J, linewidth=2.0)
ax.plot(t, R, linewidth=2.0)
ax.plot(t, D, '--', linewidth=2.0)
#ax.plot(t, C, linewidth=2.0)
labels = ['Susceptible Population', 'Exposed Individuals', 'Infectious', 'Hospitalized (Infectious)', 'Cumulative Recovered', 'Cumulative Deaths', 'Cumulative Infected']
ax.legend(labels, loc=1, fancybox=False, shadow=False)
#fig.savefig('./plots/png/' + str(count) + '.png')
fig.savefig('./plots/evd_R0=' + str(R_0) + '_incub=' + str(k) + '_fat=' + str(F) + '_tdiag=' + str(alpha) + '_isol=' + str(l) + '.pdf')
plt.close()
















