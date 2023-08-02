## Function of progesterone metabolism in cytosol by 5-alpha-reductase. 
#DOI:
def V5ared(b):
  # b is concentration of cytosolic progesterone, a is metabolite production in uM/hour. 
  km =  1.93 #$
  vmax = 34.40 #$
  a = (vmax * b)/(km + b)
  return a