## Function of monoamine transporter in vesicles.
#DOI:
# Transport of histamine from cytosol to vesicles in mast cell.
# b = cha
# c = vha
def VMATHmc(b, c):
  k = 24         
  V =  21104
  a = V*(b/(k + b)) - 5*c
  return a