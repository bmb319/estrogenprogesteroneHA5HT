## Function of mast cell travel and activation in the brain.
#DOI:
#mc_switch is the presence or lack of neuroinflammation,
# mc_start is the time of start of neuroinflammation, t is time. 
import math
def  mc_activation(t, mc_switch, mc_start_time):
  c = 1 #Amplitude of MC activation. 
  b = 12 #Strength of increase. 
  t1 = mc_start_time + math.log(999)/b #Start time is 0.01% of sigmoid function.
  if mc_switch == 0:
      f = 0
  else:
    if t < mc_start_time:
      f = 0.001
    else:
      f = c/(1 + np.exp(-b * (t - t1)))
  return f