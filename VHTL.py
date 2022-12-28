## Function of histidine transporter. 
# Histidine transport from blood to cytosol.
# UNITS in uM/h. 
# b = bht.
def VHTL(b):
  km = 1000
  vmax = 4680
  a = vmax*(b/(km + b))
  return a