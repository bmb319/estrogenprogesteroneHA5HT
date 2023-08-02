## Function of histidine decarboxylase in mast cells. 
#DOI:
# Histamine synthesis in mast cells.
# UNITS in uM/h. 
# b = gha.
def VHTDCmc(b):
# b = cht
#  c = G*
  km = 270  #homosapiens BRENDA  (other wild type =100, Komori2012)
  v = 877.50
  a =  v*(b/(km + b))
  return a