## Function of cholesterol -> pregnenolone production by P450scc enzyme. 
def VP450scc (b):
  # b is concentration of cytosolic cholesterol. a is product speed production in uM/hour. 
  km =  1 #$
  vmax = 1 #$
  a = (vmax * b)/(km + b)
  return a