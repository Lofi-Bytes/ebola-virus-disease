#***********************************************************************
# Importing everything python needs in order to be smart.
#***********************************************************************
import numpy as np                  # Numpy tells python how to do array operations.
import matplotlib.pyplot as plt     # Matplotlib tells python how to make pretty plots.
import scipy.integrate as spi       # Python needs to know how to integrate systems of differential equations.
from scipy.integrate import odeint  # Teaches Python how to solve ODE's.

def test_if():
    if t >= tint:
        betaf = beta * (1 - red/100.0)
        print beta
    else:
        betaf = beta
        print beta

t  = np.linspace(0, 500, 50000)

for i in range (0,1000):
    beta = i
    red = i + 2
    tint = 1 - 50
    test_if()


