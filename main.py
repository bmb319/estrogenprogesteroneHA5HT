from scipy.integrate._ivp.common import validate_max_step
import numpy as np
import pandas as pd
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt
import math 
from comp_model import *

#the equations model the differential terms of this model, calculate how much these variables are going to change in each iterations
#find the ks that give the diff eq to 0

#Time array.
t_factor = 1 # Time factor for graphs.
time = 10/t_factor # Time of simulation depending on t_factor.
sampling_rate = 10*t_factor #number of samples per time factor units.
time_array = np.linspace(0, time, math.floor(time * sampling_rate + 1))


## Mast cell model of neuroinflammation. 
mc_start_time = 0/t_factor #Time to start neuroinflammation effects with mast cells.
mc_switch = 0 #Switch that turns on an off all effects of mast cell presence.


#Initial conditions
z0 = [3.1968, 140.3708, 1.4717,  100, 250, 300, 0.7114, 1.3593, 1.0084, 0.7221, 1.3593, 1.0084, 150, 3,	140,	0.7205,1,1, 5.36792170e-02,  3.28039657e-03, 2.47918457e-02,  6.30086107e-04, 2.27306535e-04,  2.42721849e-04, 1] 
#Constant parameters for steady state
 
#buscar valores que tengan sentido y evidentemente sea positive


#Get solution of the differential equation.
sol = solve_ivp(comp_model, t_span = (time_array[0], time_array[-1]), t_eval = time_array, y0 = z0, method = 'RK45') #shows how the steady state values have changed over time. 



## Plot all variables.
#plt.figure(1)
#for i in range(0, 10):
#  plt.subplot(2, 5,i+1)
#  plt.plot(time_array, sol.y[i, :])
#  plt.ylabel(str(i))
#plt.show()


plt.figure(0)
plt.plot(time_array, sol.y[4, :])
plt.show()


#for i in range(0, 3):
#  plt.figure(i)
#  plt.plot(time_array, sol.y[6+i, :])
#  plt.show()

