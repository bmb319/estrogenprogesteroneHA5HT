from inhibsynHAtoHA import *
from VMATH import *
from VHNMT import *
from VHTDC import *
from VHAT import *
from VHTL import *
from HTin import *
from inhibRHAtoHA import *
from fireha import *
from mc_activation import *
from VHATmc import *
from degran_ha_mc import *
from VHTLmc import *
from VHTDCmc import *
from VMATHmc import *
from VHNMTmc import *
from activ_E2_to_ha_mc import *
from activ_E2_to_ha_syn_neuron import *
from activ_E2_to_ha_R_neuron import *
from inhib_ep_to_ha_mc import *
from inhib_cp_to_ERalpha import *
from BPin import *
from BE2in import *


import numpy as np


def comp_model_fitting(x): 
  
  dz = np.zeros(1) 
  z = np.array([100])
  q = x[0]




  b4 = 0.25
  bht0 = 100
  

  dz[0] = q*HTin(0) - VHTL(z[0])  - b4*(z[0] - bht0) 

  
  return dz