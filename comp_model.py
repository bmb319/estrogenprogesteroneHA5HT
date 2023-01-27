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
from inhib_ep_to_ha import *
import numpy as np

def comp_model(z, t):
    dz = np.zeros(len(z)) #Initialization of differential terms.
    mc_switch = 0 #Switch for the activation of mast cell neuroinflammation
    mc_start_time = 0 #Start time of mast cells activation (h)

    ## ------------------Estrogen and Progesterone in HA neuron (Soma) -------------------------------
    #Estrogen and progesterone variables.
    #z[0] Cytosolic E2. $
    #z[1] Extracellular E2. (250 pg/ml)
    #z[2] Cytosolic E1. $
    #z[3] Extracellular E1. (110 pg/ml)
    #z[4] Extracellular E3. (9 pg/ml)
    #z[5] Cytosolic E3. $
    #z[6] Cytosolic progesterone in uM. $
    #z[7]  Extracellular progesterone in uM. $
    #z[8]  Extracellular cholesterol in uM. $
    #z[9]  Cytosolic cholesterol in uM. $
    #z[10] Extracellular pregnenolone in uM. $
    #z[11]  Cytosolic pregnenolone in uM. $
    #z[12] Extracellular androstenedione in uM. $
    #z[13] Cytosolic androstenedione in uM. $
    #z[14] Extracellular testosterone in uM. $
    #z[15] Cytosolic testosterone in uM. $
    #z[16] Bound E2 to ER-alpha receptors in the nucleus of neurons in uM. $

    a1 = 1 #Constant of progesterone effect to c2. hour^-1
    a2 = 1 #Constant of histamine effects on histamine.hour^-1
    a3 = 1 #Constant of serotonin effects on serotonin. hour^-1
    a4 = 1 #Constant of blood progesterone diffusion to extracellular progesterone. hour^-1
    a5 = 1 #Constant of extracellular progesterone removal away from the synapse. hour^-1
    a6 = 1 #Constant of diffusion of progesterone from extracellular space to cytosol. hour^-1
    a7 = 1 #Constant of diffusion of blood c2 to extracellular c2. hour^-1
    a8 = 1 #Constant of diffusion of c2 from extracellular space to cytosol. hour^-1
    a9 = 1 #Constant of removal of extracellular c2 to somewhere else.
    a10 = 1 #Constant of diffusion of blood c1 to extracellular c1. hour^-1
    a11 = 1 #Constant of diffusion of c1 from extracellular space to cytosol. hour^-1
    a12 = 1 #Constant of removal of extracellular c1 to somewhere else.  hour^-1
    a13 = 1 #Constant of diffusion of c3 from blood to extracellular space. hour^-1
    a14 = 1 #Constant of diffusion of c3 from extracellular space to cytosol. hour^-1
    a15 = 1 #Constant of removal of extracellular c3 to somewhere else. hour^-1
    a16 = 1 #Constant of diffusion of cholesterol from blood to extracellular space. hour^-1
    a17 = 1 #Constant of diffusion of cholesterol from extracellular space to cytosol. hour^-1
    a18 = 1 #Constant of removal of extracellular cholesterol to somewhere else. hour^-1
    a19 = 1 #Constant of diffusion of pregnenolone from blood to extracellular space. hour^-1
    a20 = 1 #Constant of diffusion of pregnenolone from extracellular space to cytosol. hour^-1
    a21 = 1 #Constant of removal of extracellular pregnenolone to somewhere else. hour^-1
    a22 = 1 #Constant of diffusion of androstenedione from blood to extracellular space. hour^-1
    a23 = 1 #Constant of diffusion of androstenedione from extracellular space to cytosol. hour^-1
    a24 = 1 #Constant of removal of extracellular androstenedione to somewhere else. hour^-1
    a25 = 1 #Constant of diffusion of testosterone from blood to extracellular space. hour^-1
    a26 = 1 #Constant of diffusion of testosterone from extracellular space to cytosol. hour^-1
    a27 = 1 #Constant of removal of extracellular testosterone to somewhere else. hour^-1
    a28 = 0.1 #Constant of cytosol c2 binding to ER2 alpha receptor.  hour^-1
    a29 = 0.1 #Constant of bound cytosol E2 to free cytosol E2.  hour^-1
    a30 = 0.1 #Constant of free extracellular E2 to bound E2 in GPER in MC. hour^-1
    a31 = 0.1 #Constant of bound extracellular E2 in GPER to free extracellular E2 in MC. hour^-1
    a32 = 0.1 #Constant of free extracellular progesterone to bound progesterone in mPRs in MC. hour^-1
    a33 = 0.1 #Constant of bound progesterone in mPRs to free extracellular progesterone in MC. hour^-1


    # z-initial conditions and variables of differential equations.
    #Concentration values.
    bp = 1000 # Blood progesterone. uM $
    bc2 = 1000 #Blood c2.  uM. $
    bc1 = 1000 #Blood c1. uM. $
    bc3 = 10 # Blood c3. uM. $
    bch = 100 #Blood cholesterol. uM $
    bpreg = 100 #Blood pregnenolone. uM $
    bandro = 100 #Blood androstenedione in uM. $
    btest = 100 #Blood testosterone in uM. $
    basal_bound_cc2 = 1 #Basal bound cytosolic E2 to ER alpha. uM.

  
    #Equations
    dz[0] = a1*z[6] + a8 * (z[1] - z[0]) - V17betaHSD1(z[0], z[2]) + VCYP19_2(z[15]) - a28 * z[0] + a29 * z[16]

    dz[1] = a7 * (bc2 - z[1]) - a8 * (z[1] - z[0]) - a9 * z[1] 

    dz[2] = a11 * (z[3] - z[2]) + V17betaHSD1(z[0], z[2]) - VP450_3A5(z[2]) - VCOMT(z[4]) + VCYP19_1(z[13])

    dz[3] = a10 * (bc1 - z[3]) - a11 * (z[3] - z[2]) - a12 * z[3]

    dz[4] = a14 * (z[5] - z[4]) + VP450_3A5(z[2])

    dz[5] = a13 * (bc3 - z[5]) - a14 * (z[5] - z[4]) - a15 * z[5]

    dz[6] = a6 * (z[7] - z[6]) - V5ared(z[6]) + activ_E2_to_prog(z[16], basal_bound_cc2) * V3betaHSD(z[11])

    dz[7] = a4 * (bp - z[7]) - a6 * (z[7] - z[6]) - a5 * z[7]

    dz[8] = a16 * (bch - z[8]) - a17 * (z[8] - z[9]) - a18 * z[8]

    dz[9] = a17 * (z[8] - z[9]) - VP450scc(z[9])

    dz[10] = a19 * (bpreg - z[10]) - a20 * (z[10] - z[11]) - a21 * z[10]

    dz[11] = a20 * (z[10] - z[11]) + VP450scc(z[9]) - activ_E2_to_prog(z[16], basal_bound_cc2) * V3betaHSD(z[11])

    dz[12] = a22 * (bandro - z[12]) - a23 * (z[12] - z[13]) - a24 * z[12]

    dz[13] = a23 * (z[12] - z[13]) + VCYP17(z[11]) - V17betaHSD2(z[13], z[15]) - VCYP19_1(z[13])

    dz[14] = a25 * (btest - z[14]) - a26 * (z[14] - z[15]) - a27 * z[14]

    dz[15] = a26 * (z[14] - z[15]) + V17betaHSD2(z[13], z[15]) - VCYP19_2(z[15])

    dz[16] = a28 * z[0] - a29 * z[16]

    ## ------------------Histamine neuron model -------------------------------
    ## Histamine neuron variables
    #z[17] Cytosolic histamine (uM)
    #z[18] Vesicular histamine (uM)
    #z[19] Extracellular histamine (uM)
    #z[20] Blood histidine (uM)
    #z[21] Cytosolic histidine (uM)
    #z[22] Cytosolic histidine pool (uM)
    #z[23] Activated g-coupled protein H3 receptor (uM)
    #z[24] Activated T protein H3 receptor (uM)
    #z[25] Histamine bound to h3 receptor (uM)
    #z[26] Activated g-coupled protein GPER (uM)
    #z[27] Activated T protein GPER (uM)
    #z[28] E2 bound to GPER (uM)


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


    #Synapse
    dz[17] = inhibsynHAtoHA(z[23], gstar_ha_basal) * activ_E2_to_ha_syn_neuron(z[16], basal_bound_cc2) * VHTDC(z[21])  - VMATH(z[17], z[18]) -  VHNMT(z[17]) - b1*(z[17] - z[19]) + VHAT(z[19])
    dz[18] = VMATH(z[17], z[18]) - inhibRHAtoHA(z[23], gstar_ha_basal)*activ_E2_to_ha_R_neuron(z[26],gstar_E2_basal)*fireha(t)*b2*z[18]
    dz[19] = inhibRHAtoHA(z[23], gstar_ha_basal)*activ_E2_to_ha_R_neuron(z[26],gstar_E2_basal)*fireha(t)*b2*z[18] - VHAT(z[19]) + b1*(z[17] - z[19])  - b3*z[19] - mc_activation(t, mc_switch, mc_start_time) * VHATmc(z[19]) + inhibRHAtoHA(z[33], gstar_ha_basal)*degran_ha_mc(mc_activation(t, mc_switch, mc_start_time))
    dz[20] = HTin(t) - VHTL(z[20])  - b4*(z[20] - bht0) - mc_activation(t, mc_switch, mc_start_time)*VHTLmc(z[20])
    dz[21] = VHTL(z[20]) - inhibsynHAtoHA(z[23], gstar_ha_basal) * activ_E2_to_ha_syn_neuron(z[16], basal_bound_cc2) * VHTDC(z[21]) - b5*z[21] + b6*z[22]
    dz[22] = b5*z[21] - b6*z[22] - b7*z[22]
    dz[23]  = b8*z[25]**2*(g0HH - z[23]) - b9*z[24]*z[23]
    dz[24] = b10*z[23]**2*(t0HH - z[24])  - b11*z[24]
    dz[25] = b12*z[19]*(b0HH - z[25])  - b13*z[25]
    dz[26]  = b14*z[28]**2*(g0GPER - z[26]) - b15*z[27]*z[26]
    dz[27] = b16*z[26]**2*(t0GPER - z[27])  - b17*z[27]
    dz[28] =  b18*z[1]*(b0GPER - z[28])  - b19*z[28]



    ## ------------------ Mast Cells model -------------------------------
    ## Mast cell variables
    #z[29] Cytosolic histidine in mast cells (uM)
    #z[30] Cytosolic pool of histidine in mast cells (uM)
    #z[31] Cytosolic histamine in mast cells (uM)
    #z[32] Vesicular histamine in mast cells (uM)
    #z[33] Activated g-coupled protein from H3 receptors in mast cells (uM)
    #z[34] Activated T protein from H3 receptors in mast cells (uM)
    #z[35] histamine to H3 receptors in mast cells (uM)
    #z[36] E2 to GPER in mast cells (uM)
    #z[37] progesterone to mPR in mast cells (uM)

    c1 = 1 #From cHT to HTpool.
    c2 = 1 #From HTpool to cHT.
    c3 = 1 #Removal of gHT or use somewhere else.
    c4 = 100 # bound histamine to autoreceptors produce g*.
    c5 = 961.094 #T∗ facilitates the reversion of G∗ to G.
    c6 = 20 #G∗ produces T∗.
    c7 = 66.2992 #decay coefficient of T∗
    c8 = 5  #z[19] binds to autoreceptors.
    c9 = 65.6179 #z[19] dissociates from autoreceptors.
    g0Hmc = 10  #Total g-coupled protein for H3 on mast cell (uM).
    t0Hmc = 10 #Total T protein for H3 on mast cell.
    b0Hmc = 10  #Total H3 receptors on mast cell.
    basal_bound_ec2_mc = 1 #Basal bound extracellular E2 to GPER in mast cells. uM.
    basal_bound_ep_mc = 1 #Basal bound extracellular progesterone to mPRs in mast cells. uM


    dz[29] = mc_activation(t, mc_switch, mc_start_time)*VHTLmc(z[20]) - inhibsynHAtoHA(z[33], gstar_ha_basal)*VHTDCmc(z[29]) - c1*(z[29]) + c2*(z[30])
    dz[30] = c1*(z[29]) - c2*(z[30]) - c3*(z[30])
    dz[31] = inhibsynHAtoHA(z[33], gstar_ha_basal)*VHTDCmc(z[29]) - VMATHmc(z[31], z[32]) - VHNMTmc(z[31]) + mc_activation(t, mc_switch, mc_start_time) * VHATmc(z[19])
    dz[32] = VMATHmc(z[31], z[32]) - inhib_ep_to_ha(z[37], basal_bound_ep_mc)*activ_E2_to_ha_mc(z[36], basal_bound_ec2_mc)*inhibRHAtoHA(z[33], gstar_ha_basal)*degran_ha_mc(mc_activation(t, mc_switch, mc_start_time))
    dz[33] = c4*z[35]**2*(g0Hmc - z[33]) - c5*z[34]*z[33]
    dz[34] = c6*z[33]**2*(t0Hmc - z[34])  - c7*z[34]
    dz[35] = c8*z[19]*(b0Hmc - z[35]) - c9*z[35]
    dz[36] = a30 * z[1] - a31 * z[36]
    dz[37] = a32 * z[7] - a33 * z[37]

    return dz
