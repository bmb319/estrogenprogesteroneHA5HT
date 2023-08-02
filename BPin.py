## Function of replenishment of blood progesterone. 
#DOI: necessary?
#s is rate strength (h^-1), c is the constant 
def BPin(bp, bp_basal):
  s = (200)*1

  return s*(bp_basal - bp)