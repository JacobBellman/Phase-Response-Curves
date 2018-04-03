####################################        
# This script initializes Tseng's  #
#    model using the ODE class     #
####################################

#The command below is only necessary if ODE.py is up one directory
sys.path.append('..')
from ODE import *

# Declare Symbolic Variables
wc1m = sympy.Symbol('wc1m')   # wc-1 mRNA
wc2m = sympy.Symbol('wc2m')   # wc-2 mRNA
frqm = sympy.Symbol('frqm')   # frq mRNA
WC1c = sympy.Symbol('WC1c')   # WC-1 (cytoplasm)
WC2c = sympy.Symbol('WC2c')   # WC-2 (cytoplasm)
HypoFRQc = sympy.Symbol('HypoFRQc')   # Hypophosphorylated FRQ (cytoplasm)
HyperFRQc = sympy.Symbol('HyperFRQc')   # Hyperphosphorylated FRQ (cytoplasm)
HypoFRQn = sympy.Symbol('HypoFRQn')   # Hypophosphorylated FRQ (nucleus)
HyperFRQn = sympy.Symbol('HyperFRQn')   # Hyperphosphorylated FRQ (nucleus)
HypoWCCc = sympy.Symbol('HypoWCCc')   # Hypophosphorylated WCC (cytoplasm)
HyperWCCc = sympy.Symbol('HyperWCCc')   # Hyperphosphorylated WCC (cytoplasm)
HypoWCCn = sympy.Symbol('HypoWCCn')   # Hypophosphorylated WCC (nucleus)
HyperWCCn = sympy.Symbol('HyperWCCn')   # Hyperphosphorylated WCC (nucleus)
aWCC = sympy.Symbol('aWCC')   # active WCC
laWCC = sympy.Symbol('laWCC')   # light activated WCC
vvdm = sympy.Symbol('vvdm')   # vvd mRNA
VVDc = sympy.Symbol('VVDc')   # VVD (cytoplasm)
VVDn = sympy.Symbol('VVDn')   # VVD (nucleus)
WVC = sympy.Symbol('WVC')   # WCC:VVD complex
Y = sympy.Matrix([wc1m, wc2m, frqm, WC1c, WC2c, HypoFRQc, HyperFRQc,
                  HypoFRQn, HyperFRQn, HypoWCCc, HyperWCCc, HypoWCCn,
                  HyperWCCn, aWCC, laWCC, vvdm, VVDc,
                  VVDn, WVC])                         # State variables

# Model Parameters
k01 = sympy.Symbol('k01')
K01 = sympy.Symbol('K01')
H01 = sympy.Symbol('H01')
k01a = sympy.Symbol('k01a')
k02 = sympy.Symbol('k02')
k02a01 = sympy.Symbol('k02a01')
k02a02 = sympy.Symbol('k02a02')
k03 = sympy.Symbol('k03')
k03a = sympy.Symbol('k03a')
k03i = sympy.Symbol('k03i')
k04 = sympy.Symbol('k04')
k05 = sympy.Symbol('k05')
k06 = sympy.Symbol('k06')
k07 = sympy.Symbol('k07')
k08 = sympy.Symbol('k08')
k09 = sympy.Symbol('k09')
k09a = sympy.Symbol('k09a')
k10 = sympy.Symbol('k10')
k11 = sympy.Symbol('k11')
k12 = sympy.Symbol('k12')
k13 = sympy.Symbol('k13')
k14 = sympy.Symbol('k14')
k15 = sympy.Symbol('k15')
k16 = sympy.Symbol('k16')
k17 = sympy.Symbol('k17')
k18 = sympy.Symbol('k18')
k19 = sympy.Symbol('k19')
k20 = sympy.Symbol('k20')
k21 = sympy.Symbol('k21')
k22 = sympy.Symbol('k22')
k23 = sympy.Symbol('k23')
K23 = sympy.Symbol('K23')
H23 = sympy.Symbol('H23')
k24 = sympy.Symbol('k24')
k25 = sympy.Symbol('k25')
k26 = sympy.Symbol('k26')
k27 = sympy.Symbol('k27')
k28 = sympy.Symbol('k28')
k29 = sympy.Symbol('k29')
k30 = sympy.Symbol('k30')
k31 = sympy.Symbol('k31')
k32 = sympy.Symbol('k32')
k33 = sympy.Symbol('k33')
k34 = sympy.Symbol('k34')
k35 = sympy.Symbol('k35')
k36 = sympy.Symbol('k36')
k37 = sympy.Symbol('37')
k38 = sympy.Symbol('38')
k39 = sympy.Symbol('k39')

# Dictionary of sympy parameters and their values (pulsable)
p = {k01:7.3,    # k01
     K01:0.1,    # K01
     H01:4.,    # H01
     k01a:320.,   # k01a
     k02:1.19,   # k02
     k02a01:1.2,    # k02a01
     k02a02:90.,    # k02a02
     k03:1.6,    # k03
     k03a:0.03,   # k03a
     k03i:0.03,   # k03i
     k04:800.,   # k04
     k05:0.19,   # k05
     k06:0.226,  # k06
     k07:1.,     # k07
     k08:0.68,   # k08
     k09:2.,     # k09
     k09a:0.356,  # k09a
     k10:2.4,    # k10
     k11:2.5,    # k11
     k12:6.2,    # k12
     k13:0.472,  # k13
     k14:0.1,    # k14
     k15:0.3,    # k15
     k16:0.3,    # k16
     k17:0.1,    # k17
     k18:0.3,    # k18
     k19:0.29,   # k19
     k20:0.1,    # k20
     k21:0.1,    # k21
     k22:0.3,    # k22
     k23:0.6,    # k23
     K23:0.475,  # K23
     H23:12.,    # H23
     k24:0.3,    # k24
     k25:0.15,   # k25
     k26:0.,     # k26
     k27:20.,    # k27
     k28:1.8,    # k28
     k29:0.27,   # k29
     k30:0.27,   # k30
     k31:0.135,  # k31
     k32:0.085,  # k32
     k33:0.05,   # k33
     k34:0.05,   # k34
     k35:1.29,   # k35
     k36:6.,     # k36
     k37:0.24,   # k37
     k38:0.24,   # k38
     k39:0.75}   # k39

# Constant pulsing parameters
qwc1m = sympy.Symbol('qwc1m')
qwc2m = sympy.Symbol('qwc2m')
qfrqm = sympy.Symbol('qfrqm')
qWC1c = sympy.Symbol('qWC1c')
qWC2c = sympy.Symbol('qWC2c')
qHypoFRQc = sympy.Symbol('qHypoFRQc')
qHyperFRQc = sympy.Symbol('qHyperFRQc')
qHypoFRQn = sympy.Symbol('qHypoFRQn')
qHyperFRQn = sympy.Symbol('qHyperFRQn')
qHypoWCCc = sympy.Symbol('qHypoWCCc')
qHyperWCCc = sympy.Symbol('qHyperWCCc')
qHypoWCCn = sympy.Symbol('qHypoWCCn')
qHyperWCCn = sympy.Symbol('qHyperWCCn')
qaWCC = sympy.Symbol('qaWCC')
qlaWCC = sympy.Symbol('qlaWCC')
qvvdm = sympy.Symbol('qvvdm')
qVVDc = sympy.Symbol('qVVDc')
qvvdm = sympy.Symbol('qvvdm')
qVVDc = sympy.Symbol('qVVDc')

# Parameters as additive pulses           
p_add = [qwc1m, qwc2m, qfrqm, qWC1c, qWC2c,
         qHypoFRQc, qHyperFRQc, qHypoFRQn, qHyperFRQn,
         qHypoWCCc, qHyperWCCc, qHypoWCCn, qHyperWCCn,
         qaWCC, qlaWCC, qvvdm, qVVDc, qvvdm, qVVDc]

# kinetic equations
v1 = k01*aWCC**H01/(K01**H01+aWCC**H01)+k01a*laWCC
v2 = k02+k02a01*aWCC+k02a02*laWCC
v3 = k03/(1+HypoWCCn*k03i)+HypoFRQn*k03a
v4 = k04*laWCC
v5 = k05*frqm
v6 = k06*wc1m
v7 = k07*wc2m
v8 = k08*vvdm
v9 = frqm*(k09+HypoFRQc*k09a)
v10 = k10*wc1m
v11 = k11*wc2m
v12 = k12*vvdm
v13 = k13*WC1c*WC2c
v14 = k14*HypoFRQc
v15 = k15*HypoWCCc
v16 = k16*VVDc
v17 = k17*HypoFRQn
v18 = k18*HyperFRQn
v19 = k19*HyperWCCn
v20 = k20*HypoFRQc
v21 = k21*HypoFRQn
v22 = k22*HypoWCCc
v23 = k23*HypoWCCn*HypoFRQn**H23/(K23**H23+HypoFRQn**H23)
v24 = k24*HyperWCCc
v25 = k25*HypoWCCn
v26 = k26*HypoWCCn
v27 = k27*laWCC*VVDn
v28 = k28*WVC
v29 = k29*HyperFRQc
v30 = k30*HyperFRQn
v31 = k31*WC1c
v32 = k32*WC2c
v33 = k33*HyperWCCc
v34 = k34*HyperWCCn
v35 = k35*aWCC
v36 = k36*laWCC
v37 = k37*VVDc
v38 = k38*VVDn
v39 = k39*WVC

# Symbolic ODE RHS
dwc1m = v2-v10
dwc2m = v3-v11
dfrqm = v1-v9
dWC1c = v6-v31-v13
dWC2c = v7-v32-v13
dHypoFRQc = v5-v14+v17-v20
dHyperFRQc = v18+v20-v29
dHypoFRQn = v14-v17-v21
dHyperFRQn = v21-v30-v18
dHypoWCCc = v13-v15-v22+v24
dHyperWCCc = v19+v22-v24-v33
dHypoWCCn = v15-v25-v26+v28-v23
dHyperWCCn = v23-v34-v19
daWCC = v25-v35
dlaWCC = v26-v36-v27
dvvdm = v4-v12
dVVDc = v8-v37-v16
dVVDn = v16-v38-v27+v28
dWVC = v27-v39-v28

# Vector ODE RHS
dY = sympy.Matrix([dwc1m, dwc2m, dfrqm, dWC1c, dWC2c, dHypoFRQc, dHyperFRQc,
                  dHypoFRQn, dHyperFRQn, dHypoWCCc, dHyperWCCc, dHypoWCCn,
                  dHyperWCCn, daWCC, dlaWCC, dvvdm, dVVDc,
                  dVVDn, dWVC])                           
        
# Initial values on the limit cycle at frq mRNA max.
Y0 = np.array([0.5368499,  # wc1m
               0.63103175, # wc2m
               0.98507013, # frqm
               0.04024397, # WC1c
               6.09927018, # WC2c
               0.88948853, # HypoFRQc
               0.34377087, # HyperFRQc
               0.36615087, # HypoFRQn
               0.06090573, # HyperFRQn
               0.38926666, # HypoWCCc
               0.3769203,  # HyperWCCc
               0.70541941, # HypoWCCn
               0.03008159, # HyperWCCn
               0.08212058, # aWCC
               0.,         # laWCC
               0.,         # vvdm
               0.,         # VVDc
               0.,         # VVDn
               0.])        # WVC

frq_ind = 2

period = 21.62451

##### Below is example solving and plotting frq mRNA for Tseng's model #####

ModTseng = ODE('Tseng\'s Model',frq_ind,Y,p,p_add,dY,Y0,period,findQ=False)

tot_time=100
t_test = np.linspace(0,tot_time, 100*tot_time)
sol_test = ModTseng.solve(ModTseng.lim_cyc_vals, t_test)
frq_sol = sol_test[:,frq_ind]
plt.plot(t_test, frq_sol)