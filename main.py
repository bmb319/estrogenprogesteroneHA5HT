from scipy.integrate._ivp.common import validate_max_step
import numpy as np
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math 
from comp_model import *

#the equations model the differential terms of this model, calculate how much these variables are going to change in each iterations
#find the ks that give the diff eq to 0

#Time array.
t_factor = 1 # Time factor for graphs.
time = 300/t_factor # Time of simulation depending on t_factor.
sampling_rate = 10*t_factor #number of samples per time factor units.
time_array = np.linspace(0, time, math.floor(time * sampling_rate + 1))


## Mast cell model of neuroinflammation. 
mc_start_time = 0/t_factor #Time to start neuroinflammation effects with mast cells.
mc_switch = 0 #Switch that turns on an off all effects of mast cell presence.


#Initial conditions
z0 = [3.1968, 140.3708, 1.4717,  100, 250, 300, 0.7221, 1.3593, 1.0084, 0.7221, 1.3593, 1.0084, 150, 3,	140,	0.7205,	1.3539,	1.0051, 0.0593, 1, 1, 1, 1, 1, 1, 1, 1] #value at time zero SUBSTITUTE THE 1's. 
#Constant parameters for steady state
 
#buscar valores que tengan sentido y evidentemente sea positive


#Get solution of the differential equation.
x = odeint(comp_model, z0, time_array) #shows how the steady state values have changed over time


plt.figure()
plt.plot(time_array, x[:, 2])
plt.xlabel('Time (h)')
plt.ylabel('HA concentration (uM)')
plt.show()