## Function of cholesterol -> pregnenolone production by P450scc enzyme.
# https://doi.org/10.1210/jc.2011-1277
def VP450scc (b):
  # b is concentration of cytosolic cholesterol. a is product speed production in uM/hour. 
  km =  2.57 
  vmax = 1
  a = (vmax * b)/(km + b)
  return a