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
'''
def f(y, t):
    Si = y[0]
    Ei = y[1]
    Ii = y[2]
    Ji = y[3]
    Ri = y[4]
    Di = y[5]
    Ci = y[6]
    Nt = y[7]
    if t >= tint:
        betaf = beta * (1.0 - red/100.0)
    else:
        betaf = beta
    f0 = -betaf  *  Si * ((Ii + l * Ji) / Nt)                   # Susceptible individuals.
    f1 =  betaf  *  Si * ((Ii + l * Ji) / Nt) - k * Ei          # Exposed individuals.
    f2 =  k      *  Ei - (alpha  + gamma)  * Ii                 # Infectious and symptomatic individuals.
    f3 =  alpha  *  Ii - gamma_r * Ji                           # Hospitalized individuals.
    f4 =  gamma  * (1.0 - F) * Ii + gamma_r * (1.0 - F) * Ji    # Individuals removed from isolation after recovery.
    f5 =  gamma  *  F * Ii  + gamma_r * F * Ji                  # Individuals removed from isolation after disease-induced death.
    f6 =  betaf  *  Si * ((Ii + l * Ji) / Nt)                   # Cumulative second cases.
    Nt = f0 + f1 + f2 + f3 + f4                                 # Calculate the remaining population.
    return [f0, f1, f2, f3, f4, f5, f6, Nt]

'''
def f(y, t):
    Si = y[0]
    Ei = y[1]
    Ii = y[2]
    Ji = y[3]
    Ri = y[4]
    Di = y[5]
    Ci = y[6]
    Nt = y[7]
    #print 'Si_1 = ' + str(Si)
    #print 'Ei_1 = ' + str(Ei)
    #print 'Ii_1 = ' + str(Ii)
    #print 'Ji_1 = ' + str(Ji)
    #print 'Ri_1 = ' + str(Ri)
    #print 'Di_1 = ' + str(Di)
    #print 'Ci_1 = ' + str(Ci)    
    #print 'Nt_1 = ' + str(Nt)
    if t >= tint:
        betaf = beta * (1.0 - red/100.0)
    else:
        betaf = beta
    dS = -betaf  *  Si      * (Ii+l*Ji)     / Nt                      
    dE =  betaf  *  Si      * (Ii+l*Ji)     / Nt      - k       * Ei
    dI =  k      *  Ei      - (alpha+gamma) * Ii                      
    dJ =  alpha  *  Ii      - gamma_r       * Ji                      
    dR =  gamma  *  (1.0-F) * Ii            + gamma_r * (1.0-F) * Ji
    dD =  gamma  *  F       * Ii            + gamma_r * F       * Ji
    dC =  betaf  *  Si      * (Ii+l*Ji)     / Nt          
    Nt =  dS     +  dE      + dI            + dJ      + dR
    #print 'Si_2 = ' + str(Si)
    #print 'Ei_2 = ' + str(Ei)
    #print 'Ii_2 = ' + str(Ii)
    #print 'Ji_2 = ' + str(Ji)
    #print 'Ri_2 = ' + str(Ri)
    #print 'Di_2 = ' + str(Di)
    #print 'Ci_2 = ' + str(Ci)
    #print 'Nt_2 = ' + str(Nt)    
    print 'l * Ji = ' + str(l * Ji)
    print '(Ii + l * Ji) / Nt = ' + str((Ii + l * Ji) / Nt)
    print '-betaf * Si * (Ii + l * Ji) / Nt = ' + str(-betaf * Si * (Ii + l * Ji) / Nt)
    #print 'betaf = ' + str(betaf) 
    #print 'beta = ' + str(beta)
    return [dS, dE, dI, dJ, dR, dD, dC, Nt]

#***********************************************************************
# End Define Functions
#***********************************************************************

#***********************************************************************
# Define constants - Set the disease parameters.
#***********************************************************************
# True constants.
#***********************************************************************
N0      = 1000000.0             # Initial Population size.
gamma   = 1.0/5.0               # Average infectious period.
gamma_r = 1.0/7.0               # Average hospital stay.

#***********************************************************************
# Constants that need iterated.
#***********************************************************************
#R_0     = np.array([1.5, 2.0, 2.5])                                                     # Basic reprpoductive ratio of the virus.
#incub   = np.array([8.0, 10.0])                                                         # Incubation period.
#fat     = np.array([0.5, 0.6, 0.7, 0.8])                                                # Case fatality proportion.
#tdiag   = np.array([1.0, 2.0, 3.0, 4.0, 5.0])                                           # Time to diagnose.
#isol    = np.array([0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 99.9])   # Isolation effectiveness.
#tinte   = np.array([40.0, 70.0, 100.0, 500.0])                                          # Time to intervention.

R_0     = np.array([2.0 ])                                                     # Basic reprpoductive ratio of the virus.
incub   = np.array([10.0])                                                         # Incubation period.
fat     = np.array([0.5 ])                                                # Case fatality proportion.
tdiag   = np.array([2.0 ])                                           # Time to diagnose.
isol    = np.array([90.0])   # Isolation effectiveness.
tinte   = np.array([40.0])                                          # Time to intervention.
#***********************************************************************
# Set the initial conditions.
#***********************************************************************
E0 = 10.0                               # Initial number of exposed individuals.
I0 = 10.0                               # Initial number of infectious individuals in community.
J0 = 0.0                                # Initial number of hospitalized individuals.
R0 = 0.0                                # Initial number of individuals removed from isolation after recovery.
D0 = 0.0                                # Initial number of individuals removed from isolation after disease-induced death.
C0 = 0.0                                # Initial number of cumulative second cases.
S0 = N0 - (E0 + I0 + J0 + R0)           # Initial number of susceptible individuals.
y0 = [S0, E0, I0, J0, R0, D0, C0, N0]   # initial condition vector.

#***********************************************************************
# Set the minimum and maximum time in [days] to simulate.
# Set the time grid. Start at day 0 and run until day 500 and do so 
# in 5000 steps (50 steps for each day).
#***********************************************************************
t  = np.linspace(0, 500, 50000)
#t  = np.linspace(0, 1000, 100000)

#***********************************************************************
# Create a counter for naming plots sequentially (for making videos).
#***********************************************************************
count = 0

#***********************************************************************
# Begin a big-a** for loop.
#***********************************************************************
for i in range(1):
    for j in range(1):        
        for o in range(1):
            for n in range(1):
                for m in range(1):
                    for p in range(1):
                        for r in range(1):
                            #***********************************************************************
                            # Counting the reduction in transmission rate arising from control
                            # interventions or behavior changes (%).
                            #***********************************************************************
                            red = r
                            #***********************************************************************
                            # Iterate through time to intervention effect.
                            #***********************************************************************
                            #tint = tinte[p]
                            tint = 40
                            #***********************************************************************
                            # Effectiveness of isolation strategy (l=1 denotes no isolation; l=0 denotes perfect isolation).
                            # Determine l according to a given value of isol.
                            #***********************************************************************
                            #l = (100.0 - isol[m])/100.0
                            l = 0 #(100.0 - isol[m])/100.0
                            #***********************************************************************
                            # Average time from onset of symptoms to hispitalization.
                            # Determine alpha according to a given value of tdiag.
                            #***********************************************************************
                            alpha = 1.0/tdiag[n]
                            #***********************************************************************
                            # Iterate through case fatality ratios.
                            #***********************************************************************
                            F = fat[o]
                            #***********************************************************************
                            # Latency period.
                            # Determine k according to a given value of incub.
                            #***********************************************************************
                            k = 1.0/incub[j]
                            #***********************************************************************
                            # Transmission Rate.
                            # Determine beta according to a given value of R0, alpha, gamma, gamma_r, and l.
                            #***********************************************************************
                            beta = (R_0[i] * (alpha + gamma) * gamma_r) / (alpha * l + gamma_r)
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
                            C    = soln[:, 6]      # Cumulative second cases.
                            N    = soln[:, 7]      # Total population.
                            #print S[300]
                            #print E[300]
                            #print I[300]
                            #print J[300]
                            #print R[300]
                            #print D[300]
                            #print C[300]
                            #print 'N = ' + str(N[49999])
                            #***********************************************************************
                            # Update the counter for naming plots sequentially (for making videos).
                            #***********************************************************************
                            count = count + 1
                            #***********************************************************************
                            # Plot the results.
                            #***********************************************************************
                            fig = plt.figure(figsize=[12.15, 7.06])
                            ax = fig.add_subplot(111)
                            ax.set_title('Ebola Virus Outbreak [Population vs. Time]' '\n' 'R0=' + str(R_0[i]) + '  incub=' + str(incub[j]) + '  fat=' + str(fat[o]) + '  tdiag=' + str(tdiag[n]) + '  isol=' + str(isol[m]) + '  tint=' + str(tinte[p]) + '  red=' + str(red))
                            ax.set_ylabel('Population [# of People]')
                            ax.set_xlabel('Time Since Initial Outbreak [days]')
                            #ax.set_yticks(np.arange(0, 1000001, 100000))
                            #ax.semilogy(t, S, linewidth=2.0) 
                            #ax.semilogy(t, I, linewidth=2.0) 
                            #ax.semilogy(t, J, linewidth=2.0)
                            #ax.semilogy(t, R, linewidth=2.0)
                            #ax.semilogy(t, D, '--', linewidth=2.0) 
                            #ax.semilogy(t, C, linewidth=2.0)
                            #ax.semilogy(t, N, linewidth=2.0)
                            ax.plot(t, S, linewidth=2.0)
                            ax.plot(t, E, linewidth=2.0)
                            ax.plot(t, I, linewidth=2.0)
                            ax.plot(t, J, linewidth=2.0)
                            ax.plot(t, R, linewidth=2.0)
                            ax.plot(t, D, '--', linewidth=2.0)
                            ax.plot(t, C, linewidth=2.0)
                            ax.plot(t, N, linewidth=2.0)
                            labels = ['Susceptible Population', 'Exposed Individuals', 'Infectious', 'Hospitalized (Infectious)', 'Cumulative Recovered', 'Cumulative Deaths', 'Cumulative Infected', 'Total Population']
                            ax.legend(labels, loc=1, fancybox=False, shadow=False)
                            #fig.savefig('./plots/png/' + str(count) + '.png')
                            fig.savefig('./plots/evd_R0=' + str(R_0[i]) + '_incub=' + str(incub[j]) + '_fat=' + str(fat[o]) + '_tdiag=' + str(tdiag[n]) + '_isol=' + str(isol[m]) + '_tint=' + str(tinte[p]) + '_red=' + str(red) + '.pdf')
                            plt.close()
















