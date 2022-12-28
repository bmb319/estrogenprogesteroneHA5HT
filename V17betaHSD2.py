## Function of androstenedione <-> testosterone in cytosol by 17beta OH-dehydrogenase.
def V17betaHSD2(b, c):
  # b is concentration of androstenedione, c is the concentration of testosterone. 
  km1 = 1 #$
  km2 = 3.5 #uM 
  vmax_f = 1 #$
  vmax_b = 1 #$
  a =  (vmax_f * b/km1 - vmax_b*c/km2)/(1 + b/km1 + c/km1)
  return a