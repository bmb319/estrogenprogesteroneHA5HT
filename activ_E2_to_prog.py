## Function of E2 bound to ER alpha activation of progesterone synthesis. 
# a is ratio of activation, a = 1 when bound_ce2 = basal_bound_ce2.
from matplotlib import numpy


def activ_E2_to_prog(bound_ce2, basal_bound_ce2):
  c = 2 #Amplitude of activation of synthesis. 
  b = 1 # Strength of increase.
  a = c/(1 + numpy.exp(-b * (bound_ce2 - basal_bound_ce2))) 
  return a
