## Function of histamine inhibition of synthesis of histamine.
#autoreceptors and depends on levels of g-coupled proteins activated (b) by
#histamine and basal levels of g-coupled histamine (c). 
# Units in uM. 
def inhibsynHAtoHA(b, c):
  a = 1 - (1)*(b - c)
  return a