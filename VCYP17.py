## Function of pregnenolone -> androstenedione production by CYP 17 enzyme. 
def VCYP17(b):
  # b is concentration of cytosolic pregnenolone. a is product speed production in uM/hour. 
  km =  1 #$
  vmax = 1 #$
  a = (vmax * b)/(km + b)
  return a