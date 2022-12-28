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


def comp_model(z, t):

    ## --------Variables of whole model-----------------------

    #Estrogen and progesterone. 
    ce2 = z[0] #Cytosolic E2. $
    ee2 = z[1] #Extracellular E2. (250 pg/ml)
    ce1 = z[2] #Cytosolic E1. $
    ee1 = z[3] #Extracellular E1. (110 pg/ml)
    ce3 = z[4] #Extracellular E3. (9 pg/ml)
    ee3 = z[5] #Cytosolic E3. $
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
    ce2_eralpha_bound = z[16] #Bound E2 to ER-alpha receptors in the nucleus of neurons in uM. $ 


    #Histamine neuron
    cha = z[17]
    vha = z[18]
    eha = z[19]
    bht = z[20]
    cht = z[21]
    chtpool = z[22]
    gstar = z[23]
    tstar = z[24]
    bound = z[25] 

    ##Mast cell
    cht_mc = z[26] 
    chtpool_mc = z[27]
    cha_mc = z[28]
    vha_mc = z[29] 
    Gha_star_mc = z[30]
    Tha_star_mc = z[31]
    bound_ha_mc = z[32]

    mc_switch = 0
    mc_start_time = 0

    ## ------------------Estrogen and Progesterone ------------------------------- 
    k_pe = 1 #Constant of progesterone effect to E2. hour^-1
    khaha = 1 #Constant of histamine effects on histamine.hour^-1
    k5ht5ht = 1 #Constant of serotonin effects on serotonin. hour^-1
    k_bp_ep = 1 #Constant of blood progesterone diffusion to extracellular progesterone. hour^-1
    k_ep_removal = 1 #Constant of extracellular progesterone removal away from the synapse. hour^-1
    k_ep_cp = 1 #Constant of diffusion of progesterone from extracellular space to cytosol. hour^-1
    k_be2_ee2 = 1 #Constant of diffusion of blood E2 to extracellular E2. hour^-1
    k_ee2_ce2 = 1 #Constant of diffusion of E2 from extracellular space to cytosol. hour^-1
    k_ee2_removal = 1 #Constant of removal of extracellular E2 to somewhere else. 
    k_be1_ee1 = 1 #Constant of diffusion of blood E1 to extracellular E1. hour^-1
    k_ee1_ce1 = 1 #Constant of diffusion of E1 from extracellular space to cytosol. hour^-1
    k_ee1_removal = 1 #Constant of removal of extracellular E1 to somewhere else.  hour^-1
    k_be3_ee3 = 1 #Constant of diffusion of E3 from blood to extracellular space. hour^-1
    k_ee3_ce3 = 1 #Constant of diffusion of E3 from extracellular space to cytosol. hour^-1 
    k_ee3_removal = 1 #Constant of removal of extracellular E3 to somewhere else. hour^-1
    k_bch_ech = 1 #Constant of diffusion of cholesterol from blood to extracellular space. hour^-1
    k_ech_cch = 1 #Constant of diffusion of cholesterol from extracellular space to cytosol. hour^-1
    k_ch_removal = 1 #Constant of removal of extracellular cholesterol to somewhere else. hour^-1
    k_bpreg_epreg = 1 #Constant of diffusion of pregnenolone from blood to extracellular space. hour^-1
    k_epreg_cpreg = 1 #Constant of diffusion of pregnenolone from extracellular space to cytosol. hour^-1
    k_preg_removal = 1 #Constant of removal of extracellular pregnenolone to somewhere else. hour^-1
    k_bandro_eandro = 1 #Constant of diffusion of androstenedione from blood to extracellular space. hour^-1
    k_eandro_candro = 1 #Constant of diffusion of androstenedione from extracellular space to cytosol. hour^-1 
    k_andro_removal = 1 #Constant of removal of extracellular androstenedione to somewhere else. hour^-1
    k_btest_etest = 1 #Constant of diffusion of testosterone from blood to extracellular space. hour^-1
    k_etest_ctest = 1 #Constant of diffusion of testosterone from extracellular space to cytosol. hour^-1 
    k_test_removal = 1 #Constant of removal of extracellular testosterone to somewhere else. hour^-1 
    k_ce2_ce2bound = 0.1 #Constant of cytosol E2 binding to ER2 alpha receptor.  hour^-1
    k_ce2bound_ce2 = 0.1 #Constant of bound cytosol E2 to free cytosol E2.  hour^-1



    # z-initial conditions and variables of differential equations. 
    #Concentration values. 
    bp = 1000 # Blood progesterone. uM $
    be2 = 1000 #Blood E2.  uM. $ 
    be1 = 1000 #Blood E1. uM. $ 
    be3 = 10 # Blood E3. uM. $
    bch = 100 #Blood cholesterol. uM $ 
    bpreg = 100 #Blood pregnenolone. uM $
    bandro = 100 #Blood androstenedione in uM. $ 
    btest = 100 #Blood testosterone in uM. $ 
    basal_bound_ce2 = 1 #Basal bound cytosolic E2 to ER alpha. uM. 




    #Equations
    dce2 = k_pe*cp + k_ee2_ce2 * (ee2 - ce2) - V17betaHSD1(ce2, ce1) + VCYP19_2(ctest) - k_ce2_ce2bound * ce2 + k_ce2bound_ce2 * ce2_eralpha_bound

    dee2 = k_be2_ee2 * (be2 - ee2) - k_ee2_ce2 * (ee2 - ce2) - k_ee2_removal * ee2

    dce1 = k_ee1_ce1 * (ee1 - ce1) + V17betaHSD1(ce2, ce1) - VP450_3A5(ce1) - VCOMT(ce3) + VCYP19_1(candro)

    dee1 = k_be1_ee1 * (be1 - ee1) - k_ee1_ce1 * (ee1 - ce1) - k_ee1_removal * ee1

    dce3 = k_ee3_ce3 * (ee3 - ce3) + VP450_3A5(ce1) 

    dee3 = k_be3_ee3 * (be3 - ee3) - k_ee3_ce3 * (ee3 - ce3) - k_ee3_removal * ee3

    dcp = k_ep_cp * (ep - cp) - V5ared(cp) + activ_E2_to_prog(ce2_eralpha_bound, basal_bound_ce2) * V3betaHSD(cpreg)

    dep = k_bp_ep * (bp - ep) - k_ep_cp * (ep - cp) - k_ep_removal * ep

    dech = k_bch_ech * (bch - ech) - k_ech_cch * (ech - cch) - k_ch_removal * ech

    dcch = k_ech_cch * (ech - cch) - VP450scc(cch)

    depreg = k_bpreg_epreg * (bpreg - epreg) - k_epreg_cpreg * (epreg - cpreg) - k_preg_removal * epreg

    dcpreg = k_epreg_cpreg * (epreg - cpreg) + VP450scc(cch) - activ_E2_to_prog(ce2_eralpha_bound, basal_bound_ce2) * V3betaHSD(cpreg)

    deandro = k_bandro_eandro * (bandro - eandro) - k_eandro_candro * (eandro - candro) - k_andro_removal * eandro

    dcandro = k_eandro_candro * (eandro - candro) + VCYP17(cpreg) - V17betaHSD2(candro, ctest) - VCYP19_1(candro)

    detest = k_btest_etest * (btest - etest) - k_etest_ctest * (etest - ctest) - k_test_removal * etest

    dctest = k_etest_ctest * (etest - ctest) + V17betaHSD2(candro, ctest) - VCYP19_2(ctest)

    dce2_eralpha_bound = k_ce2_ce2bound * ce2 - k_ce2bound_ce2 * ce2_eralpha_bound 


    ## ------------------Histamine model -------------------------------
    b1 = 15  #HA leakage from the cytosol to the extracellular space. 
    b2 = 3.5  #HA release per action potential. 
    b4 = 0.05  #HA removal from the extracellular space
    b5 = .25  #Strength of stabilization of blood HT near 100μM. 
    b6 = 2.5 #From cHT to HTpool.
    b7 = 1 #From HTpool to cHT. 
    b8 = 1 #Other uses of HT remove HT. 
    c9 = 100 #Bound autoreceptors produceG∗. 
    c10 = 961.094 #T∗ facilitates the reversion of G∗ to G. 
    c11 = 20 #G∗ produces T∗. 
    c12 = 66.2992 #decay coefficient of T∗
    c13 = 0.5*10  #eHA binds to autoreceptors. 
    c14 = 0.5*131.23578 #eHA dissociates from autoreceptors
    g0HH = 10  #Total gstar for H3 on HA neuron
    t0HH = 10 #Total tstar for H3 on HA neuron
    b0HH = 10  #Total H3 receptors on HA neuron


    #Steady state values. 
    gstar_ha_basal =  0.7484 #Equilibrium concentration of g* histamine.
    bht0 = 100 #Steady state value of blood histidine. 



    dcha = inhibsynHAtoHA(gstar, gstar_ha_basal) * VHTDC(cht)  - VMATH(cha, vha) -  VHNMT(cha) - VDAO(cha) - b1*(cha - eha) + VHAT(eha)
    dvha = VMATH(cha, vha) - inhibRHAtoHA(gstar, gstar_ha_basal)*fireha(t)*b2*vha
    deha = inhibRHAtoHA(gstar, gstar_ha_basal)*fireha(t)*b2*vha - VHAT(eha) + b1*(cha - eha)  - b4*eha - mc_activation(t, mc_switch, mc_start_time) * VHATmc(eha) + inhibRHAtoHA(Gha_star_mc, gstar_ha_basal)*degran_ha_mc(mc_activation(t, mc_switch, mc_start_time))
    dbht = HTin(t) - VHTL(bht)  - b5*(bht - bht0) - mc_activation(t, mc_switch, mc_start_time)*VHTLmc(bht)
    dcht = VHTL(bht) - inhibsynHAtoHA(gstar, gstar_ha_basal) * VHTDC(cht) - b6*cht + b7*chtpool
    dchtpool = b6*cht - b7*chtpool - b8*chtpool
    dgstar  = c9*bound**2*(g0HH - gstar) - c10*tstar*gstar
    dtstar = c11*gstar**2*(t0HH - tstar)  - c12*tstar
    dbound = c13*eha*(b0HH - bound)  - c14*bound


    e1 = 1 #From cHT to HTpool.
    e2 = 1 #From HTpool to cHT. 
    e3 = 1 #Removal of gHT or use somewhere else. 
    e4 = 100 # Bound autoreceptors produce g*. 
    e5 = 961.094 #T∗ facilitates the reversion of G∗ to G.
    e6 = 20 #G∗ produces T∗.
    e7 = 66.2992 #decay coefficient of T∗
    e8 = 5  #eHA binds to autoreceptors. 
    e9 = 65.6179 #eHA dissociates from autoreceptors.
    g0Hmc = 10  #Total gstar for H3 on mast cell.
    t0Hmc = 10 #Total tstar for H3 on mast cell.
    b0Hmc = 10  #Total H3 receptors on mast cell.

    dcht_mc = mc_activation(t, mc_switch, mc_start_time)*VHTLmc(bht) - inhibsynHAtoHA(Gha_star_mc, gstar_ha_basal)*VHTDCmc(cht_mc) - e1*(cht_mc) + e2*(chtpool_mc)
    dchtpool_mc = e1*(cht_mc) - e2*(chtpool_mc) - e3*(chtpool_mc)
    dcha_mc = inhibsynHAtoHA(Gha_star_mc, gstar_ha_basal)*VHTDCmc(cht_mc) - VMATHmc(cha_mc, vha_mc) - VHNMTmc(cha_mc) + mc_activation(t, mc_switch, mc_start_time) * VHATmc(eha)
    dvha_mc = VMATHmc(cha_mc, vha_mc) - inhibRHAtoHA(Gha_star_mc, gstar_ha_basal)*degran_ha_mc(mc_activation(t, mc_switch, mc_start_time))
    dGha_star_mc = e4*bound_ha_mc**2*(g0Hmc - Gha_star_mc) - e5*Tha_star_mc*Gha_star_mc
    dTha_star_mc = e6*Gha_star_mc**2*(t0Hmc - Tha_star_mc)  - e7*Tha_star_mc
    dbound_ha_mc = e8*eha*(b0Hmc - bound_ha_mc) - e9*bound_ha_mc



    return [dce2, dee2, dce1, dee1, dce3, dee3, dcp, dep, dech, dcch, depreg, dcpreg, deandro, dcandro, detest, dctest, dce2_eralpha_bound, dcha, dvha, deha, dbht,  dcht, dchtpool, dgstar, dtstar, dbound,  dcht_mc, dchtpool_mc, dcha_mc, dvha_mc, dGha_star_mc, dTha_star_mc, dbound_ha_mc]

