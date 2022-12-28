## Function of E1 (estrone) -> catecholestrogens production in cytosol by Catechol-O-methyltransferase (COMT)
def VCOMT(b):
  #b is concentration of cytosolic estrone. a is product speed production in uM/hour. 
  km =  1 #$
  vmax = 1 #$
  a = (vmax * b)/(km + b)
  return a