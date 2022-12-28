##Function of androstenedione -> E1 (estrone) in cytosol by CYP19. 
def VCYP19_1(b):
  # b is concentration of cytosolic androstenedione. a is product speed production in uM/hour. 
  km =  1 #$
  vmax = 1 #$
  a = (vmax * b)/(km + b)
  return a