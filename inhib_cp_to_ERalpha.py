## Function of cytosolic progesterone decreasing ER alpha receptor binding sites. 
# cp is cytosolic progesterone while basal_cp is equilibrium cytosolic progesterone. 
# a is ratio of inhibition, a = 1 when cp = basal_cp.
from matplotlib import numpy

def inhib_cp_to_ERalpha(cp, basal_cp):
  c = 2 #Amplitude of inhibition of active sites. 
  b = 1 # Strength of decrease.
  a = c/(1 + numpy.exp(b * (cp - basal_cp))) 
  return a