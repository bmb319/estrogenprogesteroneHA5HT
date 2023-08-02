##Function of testosterone (b) -> E2 (estradiol) in cytosol by CYP19 P450 enzyme. 
#Inhibited by androstenedione (b). 
#https://doi.org/10.1016/S0021-9258(17)33694-3
def VCYP19_2(b, c):
  # b is concentration of cytosolic testosterone. a is product speed production in uM/hour. 
  km = 0.4
  ki = 0.1
  vmax = 1 #$
  k_app = km*(1 + c/ki) 
  a = (vmax * b)/(k_app + b)
  return a