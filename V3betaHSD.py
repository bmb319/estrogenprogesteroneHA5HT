## Function of pregnenolone -> progesterone production by 3-beta-HSD enzyme. 
def V3betaHSD(b):
  # b is concentration of cytosolic pregnenolone. a is product speed production in uM/hour. 
  km =  1 #$
  vmax = 1 #$
  a = (vmax * b)/(km + b)
  return a