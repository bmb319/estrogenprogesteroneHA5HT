## Function of histidine transporter. 
# Histidine transport from blood to mast cell.
#UNITS in uM/h. 
# b = bht.
def VHTLmc(b):
  km = 1000  # (lobster)(6.2-19 muM conrad05)
  vmax = 4680
  a = vmax*(b/(km + b))
  return a