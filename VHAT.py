## Function of  the putative HA transporter in HA neuron. 
# UNITS in uM/h. 
# b = eha
def  VHAT(b):
  km = 10       
  vmax = 4128.3        
  a = vmax*(b/(km + b));
  return a