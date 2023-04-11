##Function of androstenedione (b) -> E1 (estrone) in cytosol by CYP19. 
#Inhibited by testosterone (c).
#https://doi.org/10.1016/S0021-9258(17)33694-3
def VCYP19_1(b, c):
  # b is concentration of cytosolic androstenedione. a is product speed production in uM/hour. 
  km =  0.1
  ki = 0.8
  vmax = 1 #$
  k_app = km*(1 + c/ki) 
  a = (vmax * b)/(k_app + b)
  return a