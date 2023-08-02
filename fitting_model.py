from scipy.integrate._ivp.common import validate_max_step
import numpy as np
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math 
from comp_model_fitting import *
from scipy.optimize import fsolve



initial_guess = np.array([0])

sol = fsolve(comp_model_fitting, initial_guess)
print(sol[0])
