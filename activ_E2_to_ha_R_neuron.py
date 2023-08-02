## Function of E2 bounding to GPER and G* activation of histamine release in neurons.
# DOI: https://doi.org/10.3389/fendo.2020.00390
# g_ee2 is ee2 activated g* from E2 binding to GPER, basal_g_ee2 is equilibrium g*. 
from matplotlib import numpy

def activ_E2_to_ha_R_neuron(g_ee2, basal_g_ee2):
  c = 2 #Amplitude of activation of synthesis. 
  b = 1 # Strength of increase.
  a = c/(1 + numpy.exp(-b * (g_ee2 - basal_g_ee2))) 
  return a
