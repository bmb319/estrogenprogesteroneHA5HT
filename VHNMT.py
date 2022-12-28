## Function of histamine methyltransferase, 
# Histamine metabolisis in cytosol.
# UNITS in uM/h. 
# b = cha
def VHNMT(b):
  km = 4.2
  vmax = 185.5 
  a = vmax*b/(km + b)
  return a