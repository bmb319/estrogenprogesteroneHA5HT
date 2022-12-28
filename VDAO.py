## Function of diamine oxidase, 
# Histamine metabolisis in cytosol.
# UNITS in uM/h. 
# b = cha
def VDAO(b):
  km = 19
  vmax = 100
  a = vmax*b/(km + b)
  return a