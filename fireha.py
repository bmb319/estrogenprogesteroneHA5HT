## Histamine firing function. 
# Units in events/h.
import numpy as np
def fireha(t):
  f = 1
  r = 12.5
  rest = 0
  t_start = 5 + rest
  t_stop = 7 + rest
  basal = 1
  b = 2 
  time = t*3600
  if time < t_start:
      f = basal
  elif time > t_start and time < t_stop:
      f = basal + r*(1 - np.exp(-b*(time - t_start)))     
  else:
      f = basal + r*(np.exp(-b*(time - t_stop)) - np.exp(-b*time))
  return f