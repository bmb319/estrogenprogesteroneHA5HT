## Function of E2 binding to ER alpha activation of histamine synthesis in neurons.
#DOI: 10.1111/j.1471-4159.2009.06272.x
# bound ce2 is  E2 binding to ER alpha receptors, basal_bound_ce2 is equilibrium bound E2. 
from matplotlib import numpy

def activ_E2_to_ha_syn_neuron(bound_ce2, basal_bound_ce2):
  c = 2 #Amplitude of activation of synthesis. 
  b = 1 # Strength of increase.
  a = c/(1 + numpy.exp(-b * (bound_ce2 - basal_bound_ce2))) 
  return a

