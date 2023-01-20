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
from inhib_ep_to_ha import *


def comp_model(z, t):

    ## --------Variables of whole model-----------------------

    #Estrogen and progesterone. 
    cc2 = z[0] #Cytosolic c2. $
    ec2 = z[1] #Extracellular c2. (250 pg/ml)
    cc1 = z[2] #Cytosolic c1. $
    ec1 = z[3] #Extracellular c1. (110 pg/ml)
    cc3 = z[4] #Extracellular c3. (9 pg/ml)
    ec3 = z[5] #Cytosolic c3. $
    cp = z[6] #Cytosolic progesterone in uM. $
    ep = z[7] # Extracellular progesterone in uM. $
    ech = z[8] # Extracellular cholesterol in uM. $ 
    cch = z[9] # Cytosolic cholesterol in uM. $ 
    epreg = z[10] # Extracellular pregnenolone in uM. $
    cpreg = z[11] # Cytosolic pregnenolone in uM. $
    eandro = z[12] #Extracellular androstenedione in uM. $
    candro = z[13] #Cytosolic androstenedione in uM. $
    etest = z[14] #Extracellular testosterone in uM. $
    ctest = z[15] #Cytosolic testosterone in uM. $
    cc2_eralpha_bound = z[16] #Bound c2 to ER-alpha receptors in the nucleus of neurons in uM. $ 

  
    ## Histamine neuron
    cha = z[17] #Cytosolic histamine (uM)
    vha = z[18] #Vesicular histamine (uM)
    eha = z[19] #Extracellular histamine (uM) 
    bht = z[20] #Blood histidine (uM)
    cht = z[21] #Cytosolic histidine (uM)
    chtpool = z[22] #Cytosolic histidine pool (uM)
    gstar = z[23] #Activated g-coupled protein H3 receptor (uM)
    tstar = z[24] #Activated T protein H3 receptor (uM)
    bound = z[25] #Bound histamine to h3 receptor (uM)
    
    ## Mast cell
    cht_mc = z[26] #Cytosolic histidine in mast cells (uM)
    chtpool_mc = z[27] #Cytosolic pool of histidine in mast cells (uM)
    cha_mc = z[28] #Cytosolic histamine in mast cells (uM)
    vha_mc = z[29] #Vesicular histamine in mast cells (uM)
    Gha_star_mc = z[30] #Activated g-coupled protein from H3 receptors in mast cells (uM)
    Tha_star_mc = z[31] #Activated T protein from H3 receptors in mast cells (uM)
    bound_ha_mc = z[32] #Bound histamine to H3 receptors in mast cells (uM)
    bound_ec2_gper_mc = z[33] #Bound E2 to GPER in mast cells (uM)
    bound_ep_mpr_mc = z[34] #Bound progesterone to mPR in mast cells (uM)

    mc_switch = 0 #Switch for the activation of mast cell neuroinflammation
    mc_start_time = 0 #Start time of mast cells activation (h)

    ## ------------------Estrogen and Progesterone in HA neuron (Soma) ------------------------------- 
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
    a29 = 0.1 #Constant of bound cytosol c2 to free cytosol c2.  hour^-1
    a30 = 0.1 #Constant of free extracellular c2 to bound c2 in GPER in MC. hour^-1
    a31 = 0.1 #Constant of bound extracellular c2 in GPER to free extracellular c2 in MC. hour^-1
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
    basal_bound_cc2 = 1 #Basal bound cytosolic c2 to ER alpha. uM. 


    #Equations
    dcc2 = a1*cp + a8 * (ec2 - cc2) - V17betaHSD1(cc2, cc1) + VCYP19_2(ctest) - a28 * cc2 + a29 * cc2_eralpha_bound

    dec2 = a7 * (bc2 - ec2) - a8 * (ec2 - cc2) - a9 * ec2 - a30 * ec2 + a31 * bound_ec2_gper_mc

    dcc1 = a11 * (ec1 - cc1) + V17betaHSD1(cc2, cc1) - VP450_3A5(cc1) - VCOMT(cc3) + VCYP19_1(candro)

    dec1 = a10 * (bc1 - ec1) - a11 * (ec1 - cc1) - a12 * ec1

    dcc3 = a14 * (ec3 - cc3) + VP450_3A5(cc1) 

    dec3 = a13 * (bc3 - ec3) - a14 * (ec3 - cc3) - a15 * ec3

    dcp = a6 * (ep - cp) - V5ared(cp) + activ_E2_to_prog(cc2_eralpha_bound, basal_bound_cc2) * V3betaHSD(cpreg)

    dep = a4 * (bp - ep) - a6 * (ep - cp) - a5 * ep - a32 * ep + a33 * bound_ep_mpr_mc

    dech = a16 * (bch - ech) - a17 * (ech - cch) - a18 * ech

    dcch = a17 * (ech - cch) - VP450scc(cch)

    depreg = a19 * (bpreg - epreg) - a20 * (epreg - cpreg) - a21 * epreg

    dcpreg = a20 * (epreg - cpreg) + VP450scc(cch) - activ_E2_to_prog(cc2_eralpha_bound, basal_bound_cc2) * V3betaHSD(cpreg)

    deandro = a22 * (bandro - eandro) - a23 * (eandro - candro) - a24 * eandro

    dcandro = a23 * (eandro - candro) + VCYP17(cpreg) - V17betaHSD2(candro, ctest) - VCYP19_1(candro)

    detest = a25 * (btest - etest) - a26 * (etest - ctest) - a27 * etest

    dctest = a26 * (etest - ctest) + V17betaHSD2(candro, ctest) - VCYP19_2(ctest)

    dcc2_eralpha_bound = a28 * cc2 - a29 * cc2_eralpha_bound 








  
    ## ------------------Histamine neuron model -------------------------------
    b1 = 15  #HA leakage from the cytosol to the extracellular space. 
    b2 = 3.5  #HA release per action potential. 
    b3 = 0.05  #HA removal from the extracellular space
    b4 = .25  #Strength of stabilization of blood HT near 100μM. 
    b5 = 2.5 #From cHT to HTpool.
    b6 = 1 #From HTpool to cHT. 
    b7 = 1 #Other uses of HT remove HT. 
    b8 = 100 #Bound autoreceptors produceG∗. 
    b9 = 961.094 #T∗ facilitates the reversion of G∗ to G. 
    b10 = 20 #G∗ produces T∗. 
    b11 = 66.2992 #decay coefficient of T∗
    b12 = 5  #eHA binds to autoreceptors. 
    b13 = 65.61789 #eHA dissociates from autoreceptors
    g0HH = 10  #Total gstar for H3 on HA neuron
    t0HH = 10 #Total tstar for H3 on HA neuron
    b0HH = 10  #Total H3 receptors on HA neuron


    #Steady state values. 
    gstar_ha_basal =  0.7484 #Equilibrium concentration of g* histamine.
    bht0 = 100 #Steady state value of blood histidine. 


    #Synapse
    dcha = inhibsynHAtoHA(gstar, gstar_ha_basal) * VHTDC(cht)  - VMATH(cha, vha) -  VHNMT(cha) - b1*(cha - eha) + VHAT(eha)
    dvha = VMATH(cha, vha) - inhibRHAtoHA(gstar, gstar_ha_basal)*fireha(t)*b2*vha
    deha = inhibRHAtoHA(gstar, gstar_ha_basal)*fireha(t)*b2*vha - VHAT(eha) + b1*(cha - eha)  - b3*eha - mc_activation(t, mc_switch, mc_start_time) * VHATmc(eha) + inhibRHAtoHA(Gha_star_mc, gstar_ha_basal)*degran_ha_mc(mc_activation(t, mc_switch, mc_start_time))
    dbht = HTin(t) - VHTL(bht)  - b4*(bht - bht0) - mc_activation(t, mc_switch, mc_start_time)*VHTLmc(bht)
    dcht = VHTL(bht) - inhibsynHAtoHA(gstar, gstar_ha_basal) * VHTDC(cht) - b5*cht + b6*chtpool
    dchtpool = b5*cht - b6*chtpool - b7*chtpool
    dgstar  = b8*bound**2*(g0HH - gstar) - b9*tstar*gstar
    dtstar = b10*gstar**2*(t0HH - tstar)  - b11*tstar
    dbound = b12*eha*(b0HH - bound)  - b13*bound

    
  
    ## ------------------ Mast Cells model -------------------------------
    c1 = 1 #From cHT to HTpool.
    c2 = 1 #From HTpool to cHT. 
    c3 = 1 #Removal of gHT or use somewhere else. 
    c4 = 100 # Bound autoreceptors produce g*. 
    c5 = 961.094 #T∗ facilitates the reversion of G∗ to G.
    c6 = 20 #G∗ produces T∗.
    c7 = 66.2992 #decay coefficient of T∗
    c8 = 5  #eHA binds to autoreceptors. 
    c9 = 65.6179 #eHA dissociates from autoreceptors.
    g0Hmc = 10  #Total gstar for H3 on mast cell (uM).
    t0Hmc = 10 #Total tstar for H3 on mast cell.
    b0Hmc = 10  #Total H3 receptors on mast cell.
    basal_bound_ec2_mc = 1 #Basal bound extracellular c2 to GPER in mast cells. uM. 
    basal_bound_ep_mc = 1 #Basal bound extracellular progesterone to mPRs in mast cells. uM 



  
  
    dcht_mc = mc_activation(t, mc_switch, mc_start_time)*VHTLmc(bht) - inhibsynHAtoHA(Gha_star_mc, gstar_ha_basal)*VHTDCmc(cht_mc) - c1*(cht_mc) + c2*(chtpool_mc)
    dchtpool_mc = c1*(cht_mc) - c2*(chtpool_mc) - c3*(chtpool_mc)
    dcha_mc = inhibsynHAtoHA(Gha_star_mc, gstar_ha_basal)*VHTDCmc(cht_mc) - VMATHmc(cha_mc, vha_mc) - VHNMTmc(cha_mc) + mc_activation(t, mc_switch, mc_start_time) * VHATmc(eha)
    dvha_mc = VMATHmc(cha_mc, vha_mc) - inhib_ep_to_ha(bound_ep_mpr_mc, basal_bound_ep_mc)*activ_E2_to_ha_mc(bound_ec2_gper_mc, basal_bound_ec2_mc)*inhibRHAtoHA(Gha_star_mc, gstar_ha_basal)*degran_ha_mc(mc_activation(t, mc_switch, mc_start_time))
    dGha_star_mc = c4*bound_ha_mc**2*(g0Hmc - Gha_star_mc) - c5*Tha_star_mc*Gha_star_mc
    dTha_star_mc = c6*Gha_star_mc**2*(t0Hmc - Tha_star_mc)  - c7*Tha_star_mc
    dbound_ha_mc = c8*eha*(b0Hmc - bound_ha_mc) - c9*bound_ha_mc
  
    dbound_ec2_gepr_mc = a30 * ec2 - a31 * bound_ec2_gper_mc
    dbound_ep_mpr_mc = a32 * ep - a33 * bound_ep_mpr_mc

    return [dcc2, dec2, dcc1, dec1, dcc3, dec3, dcp, dep, dech, dcch, depreg, dcpreg, deandro, dcandro, detest, dctest, dcc2_eralpha_bound, dcha, dvha, deha, dbht,  dcht, dchtpool, dgstar, dtstar, dbound,  dcht_mc, dchtpool_mc, dcha_mc, dvha_mc, dGha_star_mc, dTha_star_mc, dbound_ha_mc, dbound_ec2_gepr_mc, dbound_ep_mpr_mc]

