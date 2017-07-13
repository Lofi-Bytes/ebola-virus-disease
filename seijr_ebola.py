# -*- coding: utf-8 -*-
#!/usr/bin/env python

#***********************************************************************
# Importing everything python needs in order to be smart.
#***********************************************************************
import numpy as np                  # Numpy tells python how to do array operations.
import matplotlib.pyplot as plt     # Matplotlib tells python how to make pretty plots.
import scipy.integrate as spi       # Python needs to know how to integrate systems of differential equations.
from scipy.integrate import odeint
from fractions import Fraction

#***********************************************************************
# Solving a modified SEIJR model (SEIJRDC) to simulate an Ebola virus outbreak.
#***********************************************************************
#***********************************************************************
# Define constants - Set the disease parameters.
#***********************************************************************
N       = 1000000               # Population size.
gamma   = Fraction('1/5')       # Average infectious period.
gamma_r = Fraction('1/7')       # Average hospital stay.


#R0      = 1.5                   # Basic reprpoductive ratio of the virus.
#k       = Fraction('1/10')      # Latency period.
#alpha   = Fraction('1/5')       # Average time from onset of symptoms to hispitalization.
#F       = 0.7                   # Case fatality proportion
#l       = 1.0                   # Effectiveness of isolation strategy (l=1 denotes no isolation; l=0 denotes perfect isolation)

#***********************************************************************
# Define the varying parameters as arrays:
# R0
# Incubation period
# Case fatality ratio
# Time to diagnose
# Isolation effectiveness
# Time to intervention
#***********************************************************************
'''
R0            = np.array([1.5,  2.0  , 2.5])                                         # Basic reprpoductive ratio of the virus.
incub         = np.array([8.0,  10.0 ])                                             # Incubation period of the virus.
F             = np.array([50.0, 60.0 , 70.0, 80.0])                                # Case fatality proportion.
tdiag         = np.array([1.0,  2.0  , 3.0, 4.0, 5.0])                               # Time to diagnose
isoleff       = np.array([10.0, 20.0 , 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0])  
tintervention = np.array([40.0, 70.0 , 100.0, 500.0])                              # Time to intervention.

k     = 1/incub                                                                   # Latency period.
alpha = 1/tdiag                                                                   # Average time from onset of symptoms to hispitalization.
l     = (100 - isoleff)/100                                                       # Effectiveness of isolation strategy (l=1 denotes no isolation; l=0 denotes perfect isolation)

#***********************************************************************
# Determine beta according to a given value of R0, alpha, gamma, gamma_r, and l.
#***********************************************************************
beta = R0 * (alpha + gamma) * gamma_r / (alpha * l + gamma_r)
'''

#***********************************************************************
# Test
#***********************************************************************
R0            = np.array([1.5 , 2.0  ])                                         
incub         = np.array([8.0 , 10.0 ])                                             
F             = np.array([50.0, 60.0 ])                                
tdiag         = np.array([1.0 , 2.0  ])                               
isoleff       = np.array([10.0, 20.0 ])  
tintervention = np.array([40.0, 70.0 ])                              

k     = 1/incub                                                                   
alpha = 1/tdiag                                                                   
l     = (100 - isoleff)/100                                                       
beta = np.array([1.0, 2.0])    

#***********************************************************************
# Take N 1-D sequences and return N outputs with N dimensions each, 
# such that the shape is 1 in all but one dimension and the dimension 
# with the non-unit shape value cycles through all N dimensions. 
# Use numpy.ix_ to create arrays for each dimension/parameter, and then 
# concatenate these arrays together to get the final array of 
# shape (n, d) (where n is the number of combinations, and d is the 
# number of dimensions).
#***********************************************************************
def numpy_make_grid(*args):
    # an array of ones in the overall shape, for broadcasting
    ones = np.ones([len(arg) for arg in args])
    # mesh grids of the index arrays
    midx = [(ix * ones)[None] for ix in np.ix_(*args)]
    # make into one NxD array
    idx = np.concatenate(midx).reshape((len(args), -1)).transpose()
    return idx

#***********************************************************************
# Make the n x d grid.
#***********************************************************************
grid = numpy_make_grid(R0, incub, F, tdiag, isoleff, tintervention, k, alpha, l, beta)

#***********************************************************************
# Set the initial conditions.
#***********************************************************************
S0 = N      # Initial number of susceptible individuals.
E0 = 10     # Initial number of exposed individuals.
I0 = 10     # Initial number of infectious individuals in community.
J0 = 0      # Initial number of hospitalized individuals.
R0 = 0      # Initial number of individuals removed from isolation after recovery.
D0 = 0      # Initial number of individuals removed from isolation after disease-induced death.
C0 = 0      # Initial number of cumulative second cases.

y0 = [S0, E0, I0, J0, R0, D0, C0]

#***********************************************************************
# Set the minimum and maximum time in [days] to simulate.
# Set the time grid. Start at day 0 and run until day 500 and do so 
# in 5000 steps (50 steps for each day).
#***********************************************************************
#t  = np.linspace(0, 500., 50000)
t  = np.linspace(0, 1000., 100000)

#***********************************************************************
# Set of stiffly coupled differential equations.
# Solving the system dy/dt = f(y, t), where y = [S, E, I, J, R, D, C]
#***********************************************************************
'''    
    f0 = -grid[:,9]  *    Si * (( Ii            +    grid[:,8]    * Ji )      / N  )                  
    f1 =  grid[:,9]  *    Si * (( Ii            +    grid[:,8]    * Ji )      / N  )   - grid[:,6]   * Ei    
    f2 =  grid[:,6]  *    Ei - (  grid[:,7]     +    gamma      ) * Ii               
    f3 =  grid[:,7]  *    Ii -    gamma_r       *    Ji                              
    f4 =  gamma      *  ( 1  -    grid[:,2])    *    Ii           + gamma_r   *    ( 1 - grid[:,2] ) * Ji 
    f5 =  gamma      *            grid[:,2]     *    Ii           + gamma_r   *          grid[:,2]   * Ji            
    f6 =  grid[:,9]  *    Si                    * (( Ii           + grid[:,8] * Ji )                 / N  )
'''    

def f(y, t):
    Si = y[0]
    Ei = y[1]
    Ii = y[2]
    Ji = y[3]
    Ri = y[4]
    Di = y[5]
    Ci = y[6]
    
    f0 = -beta  *  Si * ((Ii + l * Ji) / N)                # Susceptible individuals.
    f1 =  beta  *  Si * ((Ii + l * Ji) / N) - k * Ei       # Exposed individuals.
    f2 =  k     *  Ei - (alpha  + gamma)  * Ii             # Infectious and symptomatic individuals.
    f3 =  alpha *  Ii - gamma_r * Ji                       # Hospitalized individuals.
    f4 =  gamma * (1 - F) * Ii + gamma_r * (1 - F) * Ji    # Individuals removed from isolation after recovery.
    f5 =  gamma *  F * Ii  + gamma_r * F * Ji              # Individuals removed from isolation after disease-induced death.
    f6 =  beta  *  Si * ((Ii + l * Ji) / N)                # Cumulative second cases.
    
    return [f0, f1, f2, f3, f4, f5, f6]

#***********************************************************************
# Solve the Differential Equations.
#***********************************************************************
for i in range(len(grid[:,1])):
    soln = odeint(f, y0, t)
    
    S = soln[:, 0]      # Susceptible individuals.
    E = soln[:, 1]      # Exposed individuals.
    I = soln[:, 2]      # Infectious and symptomatic individuals.
    J = soln[:, 3]      # Hospitalized individuals.
    R = soln[:, 4]      # Individuals removed from isolation after recovery.
    D = soln[:, 5]      # Individuals removed from isolation after disease-induced death.
    C = soln[:, 6]      # Cumulative second cases.

    #***********************************************************************
    # Plot the results.
    #***********************************************************************
    fig1 = plt.figure(figsize=[12.15, 7.06])
    plt.plot(t, S[0], label='Population')
    plt.plot(t, I[0], label='Infectious')
    plt.plot(t, J[0], label='Hospitalized (Infectious)')
    plt.plot(t, R[0], label='Cumulative Recovered')
    plt.plot(t, D[0], label='Cumulative Deaths')
    plt.plot(t, I[0] + J[0] + C[0], label='Cumulative Infected')
    plt.xlabel('Days Since Initial Outbreak')
    plt.ylabel('Population')
    plt.title('Ebola Virus Outbreak')
    plt.legend(loc=0)
    # Save the plot and close it so that it doesn't take up all your memory.
    fig1.savefig('evd.pdf')
    plt.close()














