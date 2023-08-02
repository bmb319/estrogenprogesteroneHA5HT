## Function of E2 bounding to GPER activation of mast cells degranulation.
# DOI:  https://doi.org/10.3389/fendo.2020.00390
# a is ratio of activation, a = 1 when bound_ee2 = basal_bound_ee2.
from matplotlib import numpy

def activ_E2_to_ha_mc(bound_ee2, basal_bound_ee2):
  c = 2 #Amplitude of activation of degranulation. 
  b = 1 # Strength of increase.
  a = c/(1 + numpy.exp(-b * (bound_ee2 - basal_bound_ee2))) 
  return a

