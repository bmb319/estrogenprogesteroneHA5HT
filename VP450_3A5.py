## Function of E1 (estrone) -> E3 (estriol) in cytosol by P450 cytokines (enzymes). 
def VP450_3A5(b):
  # b is concentration of cytosolic estrone. a is product speed production in uM/hour. 
  km =  1 #$
  vmax = 1 #$
  a = (vmax * b)/(km + b)
  return a