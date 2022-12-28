## Function of histamine inhibition of histamine firing.
# Dependent on Serotonin-activated g-protein (b) and the  equilibrium  value
# of  the  activated  G-protein (c).
# Units in uM. 
def inhibRHAtoHA(b, c):
  #b = gstar
  a = 1 - (3.5)*(b - c)     
  return a