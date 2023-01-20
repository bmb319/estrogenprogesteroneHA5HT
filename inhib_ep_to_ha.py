## Function of progesterone bounding to mPRs for inhibition of mast cells degranulation.
# a is ratio of inhibition, a = 1 when bound_ee2 = basal_bound_ee2.
from matplotlib import numpy

def inhib_ep_to_ha(bound_ep, basal_bound_ep):
  c = 2 #Amplitude of inhibition of degranulation. 
  b = 1 # Strength of decrease.
  a = c/(1 + numpy.exp(b * (bound_ep - basal_bound_ep))) 
  return a