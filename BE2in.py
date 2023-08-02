## Function of replenishment of blood E2. 
#DOI: necessary?
#s is rate strength (h^-1). 
def BE2in(be2, be2_basal):
  s = (100)*1
  return s*(be2_basal - be2)