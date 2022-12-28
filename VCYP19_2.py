##Function of testosterone -> E2 (estradiol) in cytosol by CYP19. 
def VCYP19_2(b):
  # b is concentration of cytosolic testosterone. a is product speed production in uM/hour. 
  km =  1 #$
  vmax = 1 #$
  a = (vmax * b)/(km + b)
  return a