#########################################        
# This script initializes the Dovzhenok #
#       model using the ODE class       #
#########################################

#The command below is only necessary if ODE.py is up one directory
sys.path.append('..')
from ODE import *

# Sympy variables
Mf = sympy.Symbol('Mf')   # State frq mRNA
Fc = sympy.Symbol('Fc')   # State FRQ cyt
Fn = sympy.Symbol('Fn')   # State FRQ nuc
Mw = sympy.Symbol('Mw')   # State wc-1 mRNA
Wc = sympy.Symbol('Wc')   # State WC-1 cyt
Wn = sympy.Symbol('Wn')   # State WC-1 nuc
FWn = sympy.Symbol('FWn') # State FRQ:WC-1 nuc
Mc = sympy.Symbol('Mc')   # State csp-1 mRNA
C = sympy.Symbol('C')     # State CSP-1
Y = sympy.Matrix([Mf, Fc, Fn, Mw,        # State variables
                  Wc, Wn, FWn, Mc, C])

# Sympy parameters
k1 = sympy.Symbol('k1')
k2 = sympy.Symbol('k2')
k3 = sympy.Symbol('k3')
k4 = sympy.Symbol('k4')
k5 = sympy.Symbol('k5')
k6 = sympy.Symbol('k6')
k7 = sympy.Symbol('k7')
k8 = sympy.Symbol('k8')
k9 = sympy.Symbol('k9')
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
J = sympy.Symbol('J')
J1 = sympy.Symbol('J1')
J2 = sympy.Symbol('J2')
J3 = sympy.Symbol('J3')
J4 = sympy.Symbol('J4')
r = sympy.Symbol('r')
k = sympy.Symbol('k')
p_const = sympy.Symbol('p_const')

# Dictionary of sympy parameters and their values (pulsable)
p = {k1:1.8,   # k1
     k2:1.8,   # k2
     k3:0.05,  # k3
     k4:0.23,  # k4
     k5:0.27,  # k5
     k6:0.27,  # k6
     k7:0.5,   # k7
     k8:1.,    # k8
     k9:40.,   # k9
     k10:0.1,   # k10
     k11:0.05,  # k11
     k12:0.02,  # k12
     k13:50.,   # k13
     k14:1.,    # k14
     k15:5.,    # k15
     k16:0.12,  # k16
     k17:1.4,   # k17
     k18:50.,   # k18
     k19:1.4,   # k19
     J:1.25,   # J
     J1:3.,    # J1
     J2:1.,    # J2
     J3:10.,   # J3
     J4:3.,    # J4
     r:1.,     # r
     k:6.,     # k
     p_const:2.} # p_const
# Constant pulsing parameters
qMf = sympy.Symbol('qMf')
qFc = sympy.Symbol('qFc')
qFn = sympy.Symbol('qFn')
qMw = sympy.Symbol('qMw')
qWc = sympy.Symbol('qWc')
qWn = sympy.Symbol('qWn')
qFWn = sympy.Symbol('qFWn')
qMc = sympy.Symbol('qMc')
qC = sympy.Symbol('qC')

# Parameters as additive pulses           
p_add = [qMf, qFc, qFn, qMw, qWc, qWn, qFWn, qMc, qC]

# Symbolic ODE RHS
dMf = k1*(Wn**k)/(J+Wn**k)-k4*Mf
dFc = k2*Mf-(k3+k5)*Fc
dFn = k3*Fc+k14*FWn-Fn*(k6+k13*Wn) 
dMw = k7*(J1**r)/(J1**r+C**r)-k10*Mw
dWc = k8*(Fc**p_const)*Mw/((J3+Mw)*(J2+Fc**p_const))-(k9+k11)*Wc
dWn = k9*Wc-Wn*(k12+k13*Fn)+k14*FWn
dFWn = k13*Fn*Wn-(k14+k15)*FWn
dMc = k16*Wn*(J4**r)/(J4**r+C**r)-k17*Mc
dC = k18*Mc-k19*C
dY = sympy.Matrix([dMf, dFc, dFn, dMw,
                   dWc, dWn, dFWn, dMc, dC])      # Vector ODE RHS

# Initial values on the limit cycle at frq mRNA max.
Y0 = np.array([2.60537964,    # Mf0,
               8.74143086,     # Fc0
               0.0110685024,  # Fn0
               3.70361114,    # Mw0
               0.00666226484, # Wc0
               0.924361228,   # Wn0
               0.0820481701,   # FWn0
               0.05214823,  # Mc0
               1.93780501])   # C0

frq_ind = 0

period = 22.3934820284

##### Below is example solving and plotting frq mRNA for Dovzhenok's model #####

ModDovzhenok = ODE('The Dovzhenok Model',frq_ind,Y,p,p_add,dY,Y0,period,findQ=False)

tot_time=100
t_test = np.linspace(0,tot_time, 100*tot_time)
sol_test = ModDovzhenok.solve(ModDovzhenok.lim_cyc_vals, t_test)
frq_sol = sol_test[:,frq_ind]
plt.plot(t_test, frq_sol)
