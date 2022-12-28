## Function of histidine decarboxylase
# Histamine synthesis in cytosol.
# UNITS in uM/h. 
def VHTDC(b):
  # b = cht
  #  c = G*
  km = 270 
  vmax = 234
  a =  vmax*(b/(km + b))
  return a