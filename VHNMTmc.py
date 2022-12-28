## Function of histamine methyltransferase, 
# Histamine metabolisis in mast cell.
# UNITS in uM/h. 
# b = gha.
def VHNMTmc(b):
  k = 4.2  # Francis 80
  V = 21.20   
  a = (V*b/(k + b))
  return a