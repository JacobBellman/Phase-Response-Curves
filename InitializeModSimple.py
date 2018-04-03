######################################        
# This script initializes the simple #
#    model using the ODE class       #
######################################

#The command below is only necessary if ODE.py is up one directory
sys.path.append('..')
from ODE import *

# Sympy variables
A = sympy.Symbol('A')   # frq mRNA (Fm)
B = sympy.Symbol('B')   # FRQ protein (Fp)
C = sympy.Symbol('C')   # WC-1:FRQ complex (W:Fp)
D = sympy.Symbol('D')   # WC-1 protein (W)
Y = sympy.Matrix([A, B, C, D]) # State Variables

# Sympy parameters
k1 = sympy.Symbol('k1')  # Constant synthesis of WC-1 protein
k2 = sympy.Symbol('k2')  # Degradation of WC-1 protein
k3 = sympy.Symbol('k3')  # Transcription of frq mRNA
k4 = sympy.Symbol('k4')  # Degradation of frq mRNA
k5 = sympy.Symbol('k5')  # Translation of frq mRNA
k6 = sympy.Symbol('k6')  # Degradation of FRQ
k7 = sympy.Symbol('k7')  # Dissociation of WC-1:FRQ complex
k8 = sympy.Symbol('k8')  # Association of WC-1:FRQ complex
k9 = sympy.Symbol('k9')  # Degradation of WC-1:FRQ complex
K = sympy.Symbol('K')    # Threshold for frq transcription function
m = sympy.Symbol('m')    # Hill-coefficient of frq transcription function

# Dictionary of sympy parameters and their values (pulsable)
p = {k1:0.161872619,
     k2:0.01205687,
     k3:5.896782858,
     k4:0.122630429,
     k5:0.320603032,
     k6:0.040189568,
     k7:0.017518666,
     k8:0.473989872,
     k9:0.52555997,
     K:0.924,
     m:7}

# Constant pulsing parameters
qA = sympy.Symbol('qA')
qB = sympy.Symbol('qB')
qC = sympy.Symbol('qC')
qD = sympy.Symbol('qD')

# Parameters as additive pulses           
p_add = [qA, qB, qC, qD]

# Symbolic ODE RHS
dA = k3*D**m/(K**m+D**m) - k4*A     # frq mRNA
dB = k5*A - k6*B + k7*C - k8*B*D    # FRQ Protein
dC = k8*B*D - k7*C - k9*C           # WC-1:FRQ complex
dD = k1 - k2*D + k7*C - k8*B*D      # WC-1 protein
dY = sympy.Matrix([dA, dB, dC, dD]) # Vector ODE RHS

# Initial values on the limit cycle at frq mRNA max.
Y0 = np.array([0.958658,   # A0
               0.7892854,  # B0
               0.3232268,  # C0
               0.5296776]) # D0

frq_ind = 0

period = 22.

##### Below is example solving and plotting frq mRNA for the simple model #####

ModSimple = ODE('The Simple Model',frq_ind,Y,p,p_add,dY,Y0,period,findQ=False)

tot_time=100
t_test = np.linspace(0,tot_time, 100*tot_time)
sol_test = ModSimple.solve(ModSimple.lim_cyc_vals, t_test)
frq_sol = sol_test[:,frq_ind]
plt.plot(t_test, frq_sol)