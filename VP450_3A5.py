## Function of E1 (estrone) -> E3 (estriol) in cytosol by P450 cytokines (enzymes). 
# https://doi.org/10.1111/j.1365-2125.2010.03656.x
def VP450_3A5(b):
  # b is concentration of cytosolic estrone. a is product speed production in uM/hour. 
  km = 3
  vmax = 1 
  a = (vmax * b)/(km + b)
  return a