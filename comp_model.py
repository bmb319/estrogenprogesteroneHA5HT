from V5ared import *
from V17betaHSD1 import *
from VP450_3A5 import *
from VCOMT import *
from VP450scc import *
from V3betaHSD import *
from VCYP17 import *
from V17betaHSD2 import *
from VCYP19_1 import *
from VCYP19_2 import *
from activ_E2_to_prog import *
from inhibsynHAtoHA import *
from VMATH import *
from VHNMT import *
from VHTDC import *
from VHAT import *
from VHTL import *
from HTin import *
from inhibRHAtoHA import *
from fireha import *
from VDAO import *
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


import numpy as np

def comp_model(z, t):
    dz = np.zeros(len(z)) #Initialization of differential terms.
    mc_switch = 0 #Switch for the activation of mast cell neuroinflammation
    mc_start_time = 0 #Start time of mast cells activation (h)



    ## ------------------ Histamine neuron model -------------------------------
    ## Histamine neuron variables
    #z[0] Cytosolic histamine (uM)
    #z[1] Vesicular histamine (uM)
    #z[2] Extracellular histamine (uM)
    #z[3] Blood histidine (uM)
    #z[4] Cytosolic histidine (uM)
    #z[5] Cytosolic histidine pool (uM)
    #z[6] Activated g-coupled protein H3 receptor (uM)
    #z[7] Activated T protein H3 receptor (uM)
    #z[8] Histamine bound to h3 receptor (uM)
    #z[9] Activated g-coupled protein GPER (uM)
    #z[10] Activated T protein GPER (uM)
    #z[11] E2 bound to GPER (uM)


    ## Histamine neuron constants. 
    b1 = 15  #HA leakage from the cytosol to the extracellular space.
    b2 = 3.5  #HA release per action potential.
    b3 = 0.05  #HA removal from the extracellular space
    b4 = .25  #Strength of stabilization of bHT to bHT0.
    b5 = 2.5 #From cHT to HTpool.
    b6 = 1 #From HTpool to cHT.
    b7 = 1 #Other uses of HT remove HT.
    b8 = 100 #Histamine bound to autoreceptors produce G∗.
    b9 = 961.094 #T∗ facilitates the reversion of G∗ to G.
    b10 = 20 #G∗ produces T∗.
    b11 = 66.2992 #decay coefficient of T∗
    b12 = 5  #eHA binds to autoreceptors.
    b13 = 65.61789 #eHA dissociates from autoreceptors
    g0HH = 10  #Total g-coupled protein for H3 on HA neuron
    t0HH = 10 #Total T protein for H3 on HA neuron
    b0HH = 10  #Total bound H3 receptors on HA neuron
    b14 = 100 #E2 bound to GPER produce G∗ in HA neuron 
    b15 = 961.094 #T∗ in GPER facilitates the reversion of G∗ to G in HA neuron
    b16 = 20 #GPER G∗ produces T∗ in HA neuron
    b17 = 66.2992 #decay coefficient of T∗ in GPER
    b18 = (0)*5  #E2 binds to GPER receptors $
    b19 = (0)*65.61789 #E2 dissociates from GPER receptors
    g0GPER = 10  #Total g-coupled protein for GPER in HA neuron
    t0GPER = 10 #Total T protein for GPER in HA neuron
    b0GPER = 10  #Total bound E2 in GPER in HA neuron

    #Steady state values.
    gstar_ha_basal =  0.7484 #Equilibrium concentration of g* histamine in H3 receptor.
    gstar_E2_basal =  0.7484 #Equilibrium concentration of g* E2 in GPER receptor.
    bht0 = 100 #Steady state value of blood histidine.
    basal_bound_ce2 = 1


    #Synapse
    dz[0] = inhibsynHAtoHA(z[6], gstar_ha_basal) * activ_E2_to_ha_syn_neuron(z[23], basal_bound_ce2) * VHTDC(z[4])  - VMATH(z[0], z[1]) -  VHNMT(z[0]) - b1*(z[0] - z[2]) + VHAT(z[2])
    dz[1] = VMATH(z[0], z[1]) - inhibRHAtoHA(z[6], gstar_ha_basal)*activ_E2_to_ha_R_neuron(z[9], gstar_E2_basal)*fireha(t)*b2*z[1]
    dz[2] = inhibRHAtoHA(z[6], gstar_ha_basal)*activ_E2_to_ha_R_neuron(z[9], gstar_E2_basal)*fireha(t)*b2*z[1] - VHAT(z[2]) + b1*(z[0] - z[2])  - b3*z[2] - mc_activation(t, mc_switch, mc_start_time) * VHATmc(z[2]) + degran_ha_mc(mc_activation(t, mc_switch, mc_start_time))*z[15]
    dz[3] = HTin(t) - VHTL(z[3])  - b4*(z[3] - bht0) - mc_activation(t, mc_switch, mc_start_time)*VHTLmc(z[3])
    dz[4] = VHTL(z[3]) - inhibsynHAtoHA(z[6], gstar_ha_basal) * activ_E2_to_ha_syn_neuron(z[23], basal_bound_ce2) * VHTDC(z[4]) - b5*z[4] + b6*z[5]
    dz[5] = b5*z[4] - b6*z[5] - b7*z[5]
    dz[6]  = b8*z[8]**2*(g0HH - z[6]) - b9*z[7]*z[6]
    dz[7] = b10*z[6]**2*(t0HH - z[7])  - b11*z[7]
    dz[8] = b12*z[2]*(b0HH - z[8])  - b13*z[8]
    dz[9]  = b14*z[11]**2*(g0GPER - z[9]) - b15*z[10]*z[9]
    dz[10] = b16*z[9]**2*(t0GPER - z[10])  - b17*z[10]
    dz[11] =  b18*z[21]*(b0GPER - z[11])  - b19*z[11]



    ## ------------------ Mast Cells model -------------------------------
    ## Mast cell variables
    #z[12] Cytosolic histidine in mast cells (uM)
    #z[13] Cytosolic pool of histidine in mast cells (uM)
    #z[14] Cytosolic histamine in mast cells (uM)
    #z[15] Vesicular histamine in mast cells (uM)

    #z[16] bound E2 to GPER in mast cells (uM)
    #z[17] bound progesterone to mPR in mast cells (uM)

    c1 = 1 #From cHT to HTpool.
    c2 = 1 #From HTpool to cHT.
    c3 = 1 #Removal of gHT or use somewhere else.
    c4 = 100 # bound histamine to autoreceptors produce g*.
    c5 = 961.094 #T∗ facilitates the reversion of G∗ to G.
    c6 = 20 #G∗ produces T∗.
    c7 = 66.2992 #decay coefficient of T∗
    c8 = 5  #z[2] binds to autoreceptors.
    c9 = 65.6179 #z[2] dissociates from autoreceptors.
    c10 = 0.1 #Constant of free extracellular E2 to bound E2 in GPER in MC. hour^-1
    c11 = 0.1 #Constant of bound extracellular E2 in GPER to free extracellular E2 in MC. hour^-1
    c12 = 0.1 #Constant of free extracellular progesterone to bound progesterone in mPRs in MC. hour^-1
    c13 = 0.1 ##Constant of bound progesterone in mPRs to free extracellular progesterone in MC. hour^-1

    basal_bound_ec2_mc = 1 #Basal bound extracellular E2 to GPER in mast cells. uM.
    basal_bound_ep_mc = 1 #Basal bound extracellular progesterone to mPRs in mast cells. uM


    dz[12] = mc_activation(t, mc_switch, mc_start_time)*VHTLmc(z[3]) - VHTDCmc(z[12]) - c1*(z[12]) + c2*(z[13])
    dz[13] = c1*(z[12]) - c2*(z[13]) - c3*(z[13])
    dz[14] = VHTDCmc(z[12]) - VMATHmc(z[14], z[15]) - VHNMTmc(z[14]) + mc_activation(t, mc_switch, mc_start_time) * VHATmc(z[2])
    dz[15] = VMATHmc(z[14], z[15]) - inhib_ep_to_ha_mc(z[17], basal_bound_ep_mc)*activ_E2_to_ha_mc(z[16], basal_bound_ec2_mc)*degran_ha_mc(mc_activation(t, mc_switch, mc_start_time))*z[15]
    dz[16] = c10 * z[1] - c11 * z[16]
    dz[17] = c12 * z[7] - c13 * z[17]



    ## ------ Estrogen and Progesterone ----------------------------------------------
  
    #Variables
    #z[18] Blood progesterone (uM)
    #z[19] Extracellular progesterone (uM)
    #z[20] Cytosolic progesterone (uM)
    #z[20] Blood E2 (uM)
    #z[21] Extracellular E2 (uM)
    #z[22] Cytosolic E2 (uM)
    #z[23] Bound E2 to ER-alpha receptors in the nucleus of neurons (uM).
    d1 = 0.1
    d2 = 0.1
    basal_cp = 1
  

    dz[18] = 0
    dz[19] = 0
    dz[20] = 0
    dz[21] = 0
    dz[22] = 0
    dz[23] = d1 * inhib_cp_to_ERalpha(z[20], basal_cp) * z[22] - d2 * z[23]

  
  



  
    return dz
